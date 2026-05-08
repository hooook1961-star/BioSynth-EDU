import streamlit as st
import pandas as pd
import streamlit.components.v1 as components
import py3Dmol
import json
import os
import datetime
from translations import LANGUAGES
from core.chem_utils import smiles_to_3d_block, get_pubchem_data, get_chembl_data, prepare_ligand_for_docking

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

# --- ВСТАВКА: АВТОМАТИЧЕСКАЯ ОЧИСТКА ПАМЯТИ ---
if "active_smiles" not in st.session_state:
    st.session_state.active_smiles = smiles

# Если текущий smiles не совпадает с тем, что в памяти 
if st.session_state.active_smiles != smiles:
    st.session_state.prepared_pdbqt = None
    st.session_state.mol_block = None
    st.session_state.active_smiles = smiles
# ----------------------------------------------

# 4. ОСНОВНОЙ ИНТЕРФЕЙС
# --- ЗАГОЛОВОК ---
st.markdown("""
    <style>
    /* Стиль для основного заголовка */
    .responsive-title {
        font-size: 2.5rem; /* Размер для ПК по умолчанию */
        font-weight: bold;
        line-height: 1.2;
        margin-bottom: 1rem;
    }

    /* Настройки для мобильных экранов (ширина менее 768px) */
    @media (max-width: 768px) {
        .responsive-title {
            font-size: 5.5vw; /* Шрифт будет подстраиваться под ширину экрана */
        }
    }
    </style>
    """, unsafe_allow_html=True)

title_text = t.get("title_main", "BioSynth-EDU")
st.markdown(f'<p class="responsive-title">🧪 {title_text}</p>', unsafe_allow_html=True)

tab1, tab2, tab3, tab4 = st.tabs([
    t["tab_3d"], 
    t["tab_admet"], 
    t["tab_docking"], 
    t["tab_edu"]
])

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
        st.info("Место для лабораторных работ из главы 6 вашего файла")
    with itabs[3]:
        st.info("Раздел для кейсов магистрантов")