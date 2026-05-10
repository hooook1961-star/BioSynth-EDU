import streamlit as st
import pandas as pd
import openpyxl
import streamlit.components.v1 as components
from streamlit_ketcher import st_ketcher
import py3Dmol
import json
import os
import datetime
import random
from translations import LANGUAGES
from core.chem_utils import smiles_to_3d_block, get_pubchem_data, get_chembl_data, prepare_ligand_for_docking

# --- ФУНКЦИЯ ЗАГРУЗКИ ---
def get_assessment_data():
    lang = st.session_state.get('lang', 'Русский')
    
    # Маппинг колонок
    lang_map = {
        "Русский": {"q_test": "question_ru", "opt_test": "options_ru", "q_open": "Вопрос (RU)"},
        "Қазақша": {"q_test": "question_kz", "opt_test": "options_kz", "q_open": "Сұрақ (KZ)"},
        "English": {"q_test": "question_en", "opt_test": "options_en", "q_open": "Question (EN)"}
    }
    cols = lang_map.get(lang)

    try:
        # Чтение EXCEL
        df_tests = pd.read_excel('data/Тесты.xlsx')
        df_open = pd.read_excel('data/Открытые вопросы.xlsx')

        # Выборка случайных данных
        sampled_tests = df_tests.sample(n=min(10, len(df_tests))).to_dict('records')
        sampled_open = df_open.sample(n=min(3, len(df_open)))[cols['q_open']].tolist()

        return sampled_tests, sampled_open, cols
    except Exception as e:
        st.error(f"Ошибка загрузки Excel: {e}")
        return None, None, None
        
st.set_page_config(page_title="BioSynth-EDU", layout="wide")

# --- 2. ИНИЦИАЛИЗАЦИЯ SESSION STATE ---
if 'mol_block' not in st.session_state:
    st.session_state.mol_block = None
if 'lang' not in st.session_state:
    st.session_state.lang = "Русский"

# --- 3. ФУНКЦИЯ ДЛЯ ТЕКСТА ПОСОБИЯ ---
@st.cache_data
def get_chapter_text(chapter_num):
    chapters = {
        "4": "### Глава 4. Методология QSAR\nТекст из пособия...",
        "5": "### Глава 5. Кодирование структур\nОписание дескрипторов...",
        "6": "### Глава 6. Практикум\nИнструкции для Excel..."
    }
    return chapters.get(chapter_num, "Выберите главу")

# --- 4. ЛОГИКА ЯЗЫКА И ЗАГРУЗКА КАТАЛОГА ---
selected_lang = st.sidebar.selectbox(
    "🌐 Language / Тіл / Язык", 
    options=list(LANGUAGES.keys()),
    index=list(LANGUAGES.keys()).index(st.session_state.lang)
)
st.session_state.lang = selected_lang
t = LANGUAGES[st.session_state.lang]

lang_code_map = {"Русский": "ru", "Қазақша": "kz", "English": "en"}
L_CODE = lang_code_map[st.session_state.lang]

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(BASE_DIR, 'data', 'catalog.json')

@st.cache_data
def load_catalog():
    try:
        if os.path.exists(DATA_PATH):
            with open(DATA_PATH, 'r', encoding='utf-8') as f:
                return json.load(f)
        return []
    except: return []

catalog = load_catalog()

# --- 4. БОКОВАЯ ПАНЕЛЬ ---
st.sidebar.header(t["sidebar_select_mol"])

# --- ГРУППА 1: КАЗАХСТАНСКИЙ КАТАЛОГ ---
st.sidebar.subheader(t["sidebar_kaz_cat"])
kaz_options = {}
for m in catalog:
    # перевод из JSON, если его нет — стандартное имя
    display_name = m.get('name_local', {}).get(L_CODE, m.get('name', 'Unknown'))
    class_name = m.get('classification_local', {}).get(L_CODE, m.get('classification', 'Bioactiv'))
    kaz_options[f"{display_name} ({class_name})"] = m['smiles']

selected_kaz = st.sidebar.selectbox(
    t["sidebar_kaz_label"], 
    options=[t["select_placeholder"]] + list(kaz_options.keys())
)

# --- ГРУППА 2: СТАНДАРТНЫЕ ПРИМЕРЫ (ПОЛНЫЙ СПИСОК) ---
st.sidebar.subheader(t["sidebar_world_cat"])
world_examples = {
    f"Аспирин ({t['cat_analgesic']})": "CC(=O)OC1=CC=CC=C1C(=O)O",
    f"Кофеин ({t['cat_stimulant']})": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    f"Парацетамол ({t['cat_antipyretic']})": "CC(=O)NC1=CC=C(O)C=C1",
    f"Ибупрофен ({t['cat_nsaid']})": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    f"Пенициллин G ({t['cat_antibiotic']})": "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
    f"Никотин ({t['cat_alkaloid']})": "CN1CCCC1C2=CN=CC=C2",
    f"Дофамин ({t['cat_neuro']})": "C1=CC(=C(C=C1CCN)O)O"
}

selected_world = st.sidebar.selectbox(
    t["sidebar_world_label"], 
    options=[t["select_placeholder"]] + list(world_examples.keys())
)

# --- ЛОГИКА ОПРЕДЕЛЕНИЯ SMILES ---
current_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O" # Default: Aspirin

if selected_kaz != t["select_placeholder"]:
    current_smiles = kaz_options[selected_kaz]
elif selected_world != t["select_placeholder"]:
    current_smiles = world_examples[selected_world]

st.sidebar.markdown("---")
st.sidebar.header(t["sidebar_manual"])

# Итоговое определение SMILES через ввод или выбор
smiles = st.sidebar.text_input(t["sidebar_manual_label"], value=current_smiles)

# --- АВТОМАТИЧЕСКАЯ ОЧИСТКА ПАМЯТИ ---
if "active_smiles" not in st.session_state:
    st.session_state.active_smiles = smiles

# Если текущий smiles не совпадает с тем, что в памяти 
if st.session_state.active_smiles != smiles:
    st.session_state.prepared_pdbqt = None
    st.session_state.mol_block = None
    st.session_state.active_smiles = smiles

# -- АВТОМАТИЧЕСКОЕ ОПРЕДЕЛЕНИЕ current_mol ---
if 'current_mol' not in st.session_state:
    st.session_state.current_mol = None

# Обновлерие current_mol при изменении smiles
current_mol_candidate = None

# 1. Поиск в казахстанском каталоге
if selected_kaz != t["select_placeholder"]:
    current_mol_candidate = next((m for m in catalog if m['smiles'] == current_smiles), None)

# 2. Если не нашли — можно добавить поиск по другим источникам позже

if current_mol_candidate:
    st.session_state.current_mol = current_mol_candidate
elif smiles != st.session_state.get('last_smiles'):
    # Если SMILES изменился вручную — сбрасываем current_mol
    st.session_state.current_mol = None

st.session_state.last_smiles = smiles
# ----------------------------------------------

# 4. ОСНОВНОЙ ИНТЕРФЕЙС
# --- ЗАГОЛОВОК ---
st.title(f"🧪 {t.get('title_main', 'BioSynth-EDU')}")

# блок вкладок:
tab_names = [
    t.get("tab_3d", "3D"), 
    t.get("tab_admet", "ADMET"), 
    t.get("tab_docking", "Docking"), 
    t.get("tab_edu", "Education"),
    t.get("tab_project", "Project") # .get не вызывает KeyError
]
tab1, tab2, tab3, tab4, tab5 = st.tabs(tab_names)

with tab1:
    col1, col2 = st.columns([3, 1])
    
    with col1:
        c1, c2 = st.columns(2)
        if c1.button(t["btn_build_3d"], use_container_width=True):
            st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=False)
        
        if c2.button(t["btn_optimize"], use_container_width=True):
            with st.spinner(t["spinner_optimize"]):
                st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=True)

        if st.session_state.mol_block:
            # Визуализация
            view = py3Dmol.view(width=700, height=500)
            view.addModel(st.session_state.mol_block, "mol")
            view.setStyle({'stick': {'radius':0.2}, 'sphere': {'scale':0.3}})
            view.zoomTo()
            view.setBackgroundColor('#ffffff')
            components.html(view._make_html(), height=550)
            
            # 1. безопасное имя файла
            try:
                # Локализованное имя из выбора
                prefix = selected_kaz.split()[0] if selected_kaz != t["select_placeholder"] else "molecule"
            except:
                prefix = "molecule"
            
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            unique_filename = f"{prefix}_{timestamp}.sdf"

            # 2. Кнопка скачивания SDF
            st.download_button(
                label=t["download_sdf"],
                data=st.session_state.mol_block,
                file_name=unique_filename,
                mime="chemical/x-mdl-sdfile",
                help=t["download_help"]
            )
        else:
            st.info(t["info_select_mol"])

    with col2:
        st.subheader(t["header_ref"])
        data = get_pubchem_data(smiles)
        
        if data:
            st.metric(t["metric_mw"], f"{data['mw']} g/mol")
            st.metric(t["metric_logp"], data['logp'])
            st.metric(t["metric_rot_bonds"], data['rotatable_bonds'])
            
            st.divider()
            with st.spinner(t["spinner_chembl"]):
                chembl_info = get_chembl_data(data['inchikey'])
            
            if chembl_info:
                st.write(f"🧬 **ChEMBL ID:** `{chembl_info['chembl_id']}`")
                
                phase = chembl_info['max_phase']
                status_color = "green" if phase == 4 else "orange"
                
                # Локализация статуса
                status_text = f"Phase {phase}"
                if phase == 4:
                    status_text += f" ({t['status_approved']})"
                
                st.markdown(f"**{t['status_label']}** <span style='color:{status_color}'>{status_text}</span>", unsafe_allow_html=True)
                
                with st.expander(t["mechanism_label"]):
                    for m in chembl_info['mechanisms']:
                        st.write(f"• {m}")
            else:
                st.caption(t["no_chembl"])

            st.divider()
            st.write(t["ext_links"])
            
            pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/#query={data['inchikey']}"
            st.link_button(t["btn_pubchem"], pubchem_url, use_container_width=True)
            
            chembl_url = f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={data['inchikey']}"
            st.link_button(t["btn_chembl"], chembl_url, use_container_width=True)
            
            chembl_sim_url = f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={smiles}&search_type=similarity&similarity=70"
            st.link_button(t["btn_similarity"], chembl_sim_url, use_container_width=True, type="primary")

        else:
            st.warning(t["warn_no_pubchem"])

with tab2:
    st.header(t["admet_header"])
    
    # Блок с инструкциями
    st.markdown(t["admet_instructions"])
    
    # Кнопки со ссылками
    col_btn1, col_btn2 = st.columns(2)
    with col_btn1:
        st.link_button(t["btn_open_admetlab"], "https://admetlab3.scbdd.com/", use_container_width=True)
    with col_btn2:
        st.link_button(t["btn_open_swissadme"], "http://www.swissadme.ch/", use_container_width=True)
    
    st.divider()
    
    # Загрузка файла
    uploaded_file = st.file_uploader(t["uploader_label"], type="csv")
    
    if uploaded_file:
        try:
            df = pd.read_csv(uploaded_file)
            
            with st.expander(t["expander_raw"]):
                st.dataframe(df)
                
            st.subheader(t["subheader_interp"])
            cols = df.columns.tolist()
            c1, c2, c3 = st.columns(3)

            # (Функция safe_float)

            # 1. Липофильность (LogP)
            logp_col = next((c for c in cols if 'logp' in c.lower()), None)
            if logp_col:
                val = safe_float(df[logp_col].iloc[0])
                with c1:
                    st.metric(t["metric_lipophilicity"], f"{val:.2f}")
                    if -1 < val < 5: 
                        st.success(t["status_optimal"])
                    else:
                        st.warning(t["status_extreme"])

            # 2. Гематоэнцефалический барьер (BBB)
            bbb_col = next((c for c in cols if 'bbb' in c.lower()), None)
            if bbb_col:
                raw_val = df[bbb_col].iloc[0]
                val = safe_float(raw_val)
                with c2:
                    st.metric(t["metric_bbb"], str(raw_val))
                    if val > 0.5 or str(raw_val).lower() == 'yes':
                        st.warning(t["status_bbb_yes"])
                    else:
                        st.success(t["status_bbb_no"])

            # 3. Токсичность (hERG)
            herg_col = next((c for c in cols if 'herg' in c.lower() or 'tox' in c.lower()), None)
            if herg_col:
                raw_val = df[herg_col].iloc[0]
                val = safe_float(raw_val)
                with c3:
                    st.metric(t["metric_tox"], str(raw_val))
                    if val > 0.5 or str(raw_val).lower() == 'yes':
                        st.error(t["status_tox_high"])
                    else:
                        st.success(t["status_tox_low"])

            st.divider()
            
            # --- Анализ правил Липинского ---
            st.subheader(t["header_lipinski"])
            
            # ---Логика подсчета violations ---
            
            if violations == 0:
                st.balloons()
                st.success(t["lipinski_success"])
            else:
                st.warning(f"{t['lipinski_warn']}{violations}.")
                st.info(t["lipinski_info"])

        except Exception as e:
            st.error(f"{t['error_interp']}{e}")
            
# --- Вкладка Докинг ---            
with tab3:
    st.header(t["docking_header"])
    
    if st.session_state.get('mol_block') or smiles:
        st.success(t["docking_mol_ready"])
        
        col_prep1, col_prep2 = st.columns(2)
        
        with col_prep1:
            st.markdown(t["docking_checklist_title"])
            st.markdown(t["docking_checklist_items"])
            
            if st.button(t["btn_run_prep"], use_container_width=True):
                with st.spinner(t["spinner_meeko"]):
                    # Подготовка данных
                    pdbqt_data = prepare_ligand_for_docking(smiles)
                    if pdbqt_data:
                        st.session_state.prepared_pdbqt = pdbqt_data
                        # Обновление 3D-сетки для корректного отображения
                        st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=True)
                        st.balloons()
                    else:
                        st.error(t["docking_error"])

        with col_prep2:
            st.info(t["docking_note_student"])
            
            # Отображение результатов при наличии в сессии
            if st.session_state.get('prepared_pdbqt'):
                st.download_button(
                    label=t["btn_download_pdbqt"],
                    data=st.session_state.prepared_pdbqt,
                    file_name="ligand.pdbqt",
                    mime="text/plain",
                    use_container_width=True
                )
                
                st.write(t["docking_view_label"])
                
                # Визуализация SDF-блока (наиболее точного для связей)
                try:
                    view = py3Dmol.view(width=400, height=400)
                    view.addModel(st.session_state.mol_block, "sdf")
                    view.setStyle({'stick': {'color': 'spectrum', 'radius': 0.15}, 'sphere': {'scale': 0.25}})
                    view.zoomTo()
                    view.setBackgroundColor('#ffffff')
                    
                    import streamlit.components.v1 as components
                    components.html(view._make_html(), height=410)
                    st.caption(t["docking_caption"])
                except Exception as e:
                    st.error(f"Render error: {e}")
            
        st.divider()
        st.subheader(t["docking_next_steps_header"])
        st.write(t["docking_next_steps_text"])
        
    else:
        st.warning(t["docking_warn_no_3d"])

# --- Вкладка Обучение ---
with tab4:
    # --- ОСТАВЛЯЕМ ВАШ БЛОК С КАЗАХСТАНСКИМИ ПРЕПАРАТАМИ (БЕЗ ИЗМЕНЕНИЙ) ---
    current_mol = next((m for m in catalog if m['smiles'] == smiles), None)
    if current_mol:
        with st.expander("🇰🇿 Сведения о казахстанской разработке", expanded=True):
            st.markdown(f"### Препарат: {current_mol['name']}")
            st.write(f"**Авторы:** {', '.join(current_mol.get('authors', []))}")
            st.info(f"**Краткое описание:** {current_mol.get('description', '')}")
        st.markdown("---") 

    # --- НОВИНКА: КНОПКИ-ВКЛАДКИ ДЛЯ РАЗДЕЛЕНИЯ КОНТЕНТА ---
    # Это создаст 4 кнопки сверху, чтобы не было "каши" на одной странице
    itabs = st.tabs(["📖 Пособие", "🎥 Видео", "🧪 Лабораторные", "🔬 Магистрантам"])

    # РАЗДЕЛ 1: ПОСОБИЕ
    with itabs[0]:
        c1, c2 = st.columns([2, 1]) # Делим экран: слева текст, справа бот
        with c1:
            ch = st.selectbox("Выберите главу:", ["4", "5", "6"], key="ch_sel")
            st.markdown(get_chapter_text(ch)) # Выводим текст из функции сверху
        with c2:
            st.subheader("ИИ-Тьютор")
            st.text_input("Вопрос по пособию:", key="bot_q")

    # РАЗДЕЛ 2: ВИДЕО (то, что у вас уже было)
    with itabs[1]:
        st.video("https://www.youtube.com/watch?v=dQw4w9WgXcQ")
        st.radio("Тест: Как изменится LogP?", ["Увеличится", "Уменьшится"], key="v_test")

    # РАЗДЕЛ 3 и 4: ПОКА ПУСТЫЕ (ДЛЯ ЛАБ И КЕЙСОВ)
    with itabs[2]:
        st.info("Лабораторные работы")
    with itabs[3]:
        st.info("Кейсы")
        
# --- ВКЛАДКА 5: ИССЛЕДОВАТЕЛЬСКИЙ ПРОЕКТ ---
with tab5:
    st.header(t.get("tab_project", "🚀 Исследовательский проект"))

    # -- ОПИСАНИЕ ПРОЕКТА ---
    st.markdown(t.get(
        "project_desc",
        "**Исследовательский проект по органической химии для студентов и учащихся**"
    ))

    st.markdown(t.get(
        "project_start_hint",
        "👈 Для начала работы выберите молекулу в боковой панели\n(раздел «База kz» или поле ручного ввода SMILES)"
    ))

    st.divider()

    mol_data = st.session_state.get('current_mol')
    project_smiles = mol_data.get('smiles', smiles) if mol_data else smiles

    if project_smiles and project_smiles.strip():

        # 1. Название и SMILES
        if mol_data:
            st.subheader(
                f"🔬 {mol_data.get('name', t.get('selected_mol', 'Выбранное соединение'))}"
            )

            if mol_data.get('name_local'):
                st.caption(
                    mol_data.get('name_local', {}).get(L_CODE, '')
                )

        else:
            st.subheader(t.get("input_struct", "🔬 Введённая структура"))

        st.info(
            f"**{t.get('current_smiles', 'Текущий SMILES')}:** `{project_smiles}`"
        )

        # 2. Информация из каталога
        if mol_data:
            st.markdown(f"### {t.get('kz_info', '🇰🇿 Сведения о казахстанской разработке')}")

            col_info1, col_info2 = st.columns([2, 1])

            with col_info1:
                st.write(
                    f"**{t.get('authors', 'Авторы')}:** "
                    f"{', '.join(mol_data.get('authors', ['—']))}"
                )

                st.info(
                    f"**{t.get('description', 'Описание')}:** "
                    f"{mol_data.get('description', t.get('info_missing', 'Информация отсутствует'))}"
                )

            with col_info2:
                if mol_data.get('classification'):
                    st.metric(
                        t.get('classification', 'Классификация'),
                        mol_data.get('classification', '—')
                    )

                if mol_data.get('year'):
                    st.metric(
                        t.get('year', 'Год'),
                        mol_data.get('year', '—')
                    )

            st.divider()

        # 3. Задание
        st.subheader(
            t.get("task_title", "📝 Задание на исследовательский проект")
        )

        with st.expander(
            t.get("task_full", "Открыть полное задание"),
            expanded=True
        ):
            st.markdown(t.get(
                "task_step_1",
                "1. Проведите прогноз спектра биологической активности с помощью **PASS Online**."
            ))

            st.markdown(t.get(
                "task_step_2",
                "2. Оцените фармакологические свойства (во вкладке **ADMET**) и проверьте выполнение правила Липинского."
            ))

            st.markdown(t.get(
                "task_step_3",
                "3. В редакторе ниже **измените структуру молекулы** и проанализируйте изменения."
            ))

            st.markdown(t.get(
                "task_step_4",
                "4. Сформулируйте выводы и рекомендации для доклада."
            ))

        st.divider()

        # 4. Редактор структуры
        st.subheader(
            t.get("editor_title", "🧪 Редактор структуры молекулы")
        )

        st.markdown(
            f"**{t.get('editor_task', 'Задание: Измените структуру ниже и нажмите «Применить изменения».')}**"
        )

        editor_key = (
            f"project_ketcher_"
            f"{hash(project_smiles) % 100000}_"
            f"{st.session_state.lang}"
        )

        edited = st_ketcher(
            project_smiles,
            key=editor_key
        )

        if edited and edited != project_smiles:

            st.success(
                t.get("structure_changed", "✅ Структура изменена в редакторе!")
            )

            st.info(
                f"**{t.get('new_smiles', 'Новый SMILES')}:**\n`{edited}`"
            )

            col_btn1, col_btn2 = st.columns(2)

            with col_btn1:
                if st.button(
                    t.get(
                        "apply_btn",
                        "🔄 Применить изменения и обновить проект"
                    ),
                    use_container_width=True,
                    type="primary"
                ):

                    if mol_data:
                        st.session_state.current_mol['smiles'] = edited

                    st.session_state.active_smiles = edited

                    st.success(
                        t.get("change_applied", "✅ Изменения применены!")
                    )

                    st.rerun()

            with col_btn2:
                st.download_button(
                    label=t.get(
                        "download_btn",
                        "💾 Скачать изменённый SMILES"
                    ),
                    data=edited,
                    file_name=(
                        f"{t.get('download_filename', 'modified_molecule')}_"
                        f"{datetime.datetime.now().strftime('%H%M')}.smi"
                    ),
                    mime="text/plain",
                    use_container_width=True
                )

        else:
            st.caption(
                t.get(
                    "change_editor_hint",
                    "👆 Измените молекулу в редакторе выше"
                )
            )

        st.divider()

        # --- 5. ОБУЧАЮЩИЕ МАТЕРИАЛЫ И СЕРВИСЫ ---
        st.subheader("📚 " + t.get(
            "tools_header",
            "Обучение и инструменты"
        ))

        # Лекции
        st.markdown(
            f"**📺 {t.get('lectures_header', 'Видео-лекции по теме')}:**"
        )

        l_col1, l_col2, l_col3 = st.columns(3)

        with l_col1:
            st.link_button(
                t.get("pass_lecture", "🎥 Лекция: PASS"),
                "https://youtube.com/...",
                use_container_width=True,
                type="secondary"
            )

        with l_col2:
            st.link_button(
                t.get("adme_lecture", "🎥 Лекция: ADME"),
                "https://youtube.com/...",
                use_container_width=True,
                type="secondary"
            )

        with l_col3:
            st.link_button(
                t.get("pubchem_lecture", "🎥 Лекция: PubChem"),
                "https://youtube.com/...",
                use_container_width=True,
                type="secondary"
            )

        # Инструменты
        st.markdown(
            f"**🛠 {t.get('tools_header', 'Инструменты анализа')}:**"
        )

        s_col1, s_col2, s_col3 = st.columns(3)

        with s_col1:
            st.link_button(
                t.get("pass_online", "🌐 PASS Online"),
                "http://www.way2drug.com/passonline/",
                use_container_width=True,
                type="primary"
            )

        with s_col2:
            st.link_button(
                t.get("swissadme", "🧪 SwissADME"),
                "http://www.swissadme.ch/",
                use_container_width=True,
                type="primary"
            )

        with s_col3:
            st.link_button(
                t.get("pubchem_search", "📊 PubChem Search"),
                f"https://pubchem.ncbi.nlm.nih.gov/#query={project_smiles}",
                use_container_width=True,
                type="primary"
            )

        # --- ТЕСТ ---
        st.divider()

        if st.button(
            t.get(
                "quiz_btn",
                "📝 Пройти тест и получить вопросы к защите"
            ),
            use_container_width=True,
            type="primary"
        ):

            tests, open_qs, cols = get_assessment_data()

            if tests:

                dialog_title = t.get(
                    "quiz_title",
                    "Интеллектуальный тренажер BioSynth-EDU"
                )
                
                @st.dialog(dialog_title, width="large")
                def run_quiz_dialog(t_data, o_qs, c_map):

                    st.write(
                        "### " + t.get(
                            "quiz_part_1",
                            "Часть 1: Тестирование"
                        )
                    )

                    with st.form("quiz_form"):

                        user_answers = []

                        for i, item in enumerate(t_data):

                            col_q = c_map['q_test']
                            col_opt = c_map['opt_test']

                            q_text = item[col_q]

                            raw_options_str = str(item[col_opt])

                            raw_options = [
                                opt.strip()
                                for opt in raw_options_str.split(';')
                            ]

                            correct_answer = raw_options[0]

                            state_key = f"shuffled_{i}_{col_opt}"
                            radio_key = f"quiz_radio_{i}_{col_opt}"

                            if state_key not in st.session_state:
                                st.session_state[state_key] = random.sample(
                                    raw_options,
                                    len(raw_options)
                                )

                            st.write(f"**{i+1}. {q_text}**")

                            ans = st.radio(
                                f"{t.get('question_label', 'Вопрос')} {i}",
                                options=st.session_state[state_key],
                                key=radio_key,
                                index=None,
                                label_visibility="collapsed"
                            )

                            user_answers.append(
                                (ans, correct_answer)
                            )

                            st.divider()

                        submit_quiz = st.form_submit_button(
                            t.get(
                                "check_result",
                                "Проверить результат"
                            ),
                            use_container_width=True
                        )

                    if submit_quiz:

                        score = sum(
                            1
                            for ans, correct in user_answers
                            if ans == correct
                        )

                        st.success(
                            f"{t.get('quiz_score', 'Ваш результат')}: "
                            f"{score} / {len(t_data)}"
                        )

                        st.write(
                            "### " + t.get(
                                "quiz_part_2",
                                "Часть 2: Вопросы для подготовки к докладу"
                            )
                        )

                        for q in o_qs:
                            st.info(q)

                        if st.button(
                            t.get("close", "Закрыть")
                        ):

                            for i in range(len(t_data)):

                                state_key = (
                                    f"shuffled_{i}_{c_map['opt_test']}"
                                )

                                if state_key in st.session_state:
                                    del st.session_state[state_key]

                            st.rerun()

                run_quiz_dialog(tests, open_qs, cols)

    else:

        st.warning(
            t.get(
                "mol_not_selected",
                "⚠️ Молекула не выбрана"
            )
        )

        st.info(
            f"""
**{t.get('how_start', 'Как начать:')}**

{t.get('start_step_1', '• Выберите соединение из Казахстанского каталога в боковой панели')}

{t.get('start_step_2', '• Или введите SMILES вручную в боковой панели')}
"""
        )
        
# --- МОДУЛЬ ИИ-ТЬЮТОРА (ВСТАВИТЬ В КОНЕЦ ФАЙЛА) ---

def load_tutor_knowledge():
    """Функция для поиска и загрузки JSON-базы из корня проекта"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Выходим на уровень выше (в корень), так как data лежит там
    root_dir = os.path.dirname(current_dir)
    file_path = os.path.join(root_dir, 'data', 'bot_knowledge_new.json')
    
    # Альтернативный путь (если запускаете локально и data внутри app)
    alt_path = os.path.join(current_dir, 'data', 'bot_knowledge_new.json')
    
    final_path = file_path if os.path.exists(file_path) else alt_path

    if os.path.exists(final_path):
        try:
            with open(final_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            st.error(f"Ошибка чтения базы: {e}")
            return None
    return None

@st.dialog("🤖 Ассистент BioSynth-EDU")
def tutor_dialog():
    """Окно чата с Тьютором"""
    kb = load_tutor_knowledge()
    
    if kb is None:
        st.error("Критическая ошибка: База данных 'bot_knowledge_new.json' не найдена в папке /data/")
        st.info("Проверьте, что папка data находится в корне вашего репозитория на GitHub.")
        return

    st.write("Привет! Я помогу с навигацией по кнопкам и объясню теорию.")
    
    if "tutor_chat_history" not in st.session_state:
        st.session_state.tutor_chat_history = []

    chat_container = st.container(height=400)
    
    for msg in st.session_state.tutor_chat_history:
        chat_container.chat_message(msg["role"]).write(msg["content"])

    if prompt := st.chat_input("Например: 'Где вкладка ADMET?'"):
        st.session_state.tutor_chat_history.append({"role": "user", "content": prompt})
        chat_container.chat_message("user").write(prompt)
        
        prompt_l = prompt.lower()
        nav = kb.get('interface_navigation_and_logic', {})
        
        # 1. Логика навигации
        if any(word in prompt_l for word in ["где", "кнопка", "найти", "как"]):
            if "3d" in prompt_l:
                tab = nav.get('Tabs_Functional', {}).get('Tab_3D', {})
                res = f"📍 **Навигация:** Вкладка **'{tab.get('title')}'**. Кнопки: `{', '.join(tab.get('buttons', []))}`."
            elif "admet" in prompt_l:
                tab = nav.get('Tabs_Functional', {}).get('Tab_ADMET', {})
                res = f"📍 **Навигация:** Вкладка **'{tab.get('title')}'**. Загрузка через `{tab.get('buttons', [0])}`."
            else:
                res = "Я знаю путь к вкладкам 3D, ADMET, PASS и Докинга. Что именно показать?"
        
        # 2. Логика теории
        elif "pass" in prompt_l or "pa" in prompt_l:
            val = kb.get('PASS_Online_Full_Knowledge_Base', {}).get('Validation_and_Robustness', {})
            res = f"🧠 **Теория:** PASS имеет точность {val.get('Accuracy', '95%')} и устойчив к потере данных."
        else:
            res = "Я готов помочь! Спросите о кнопках интерфейса или теории из курса."

        chat_container.chat_message("assistant").write(res)
        st.session_state.tutor_chat_history.append({"role": "assistant", "content": res})

# --- КНОПКА В УГЛУ ЭКРАНА ---
st.markdown("""
    <style>
    .floating-tutor {
        position: fixed;
        bottom: 30px;
        right: 30px;
        z-index: 1000;
    }
    .floating-tutor button {
        width: 70px !important;
        height: 70px !important;
        border-radius: 50% !important;
        background-color: #2e7d32 !important;
        color: white !important;
        font-size: 30px !important;
        box-shadow: 0px 4px 15px rgba(0,0,0,0.4) !important;
    }
    </style>
    """, unsafe_allow_html=True)

with st.container():
    st.markdown('<div class="floating-tutor">', unsafe_allow_html=True)
    if st.button("🤖", key="tutor_trigger_btn"):
        tutor_dialog()
    st.markdown('</div>', unsafe_allow_html=True)
