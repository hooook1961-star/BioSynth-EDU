import streamlit as st
import pandas as pd
import openpyxl
import numpy as np
import streamlit.components.v1 as components
from streamlit_ketcher import st_ketcher
import py3Dmol
import json
import os
import datetime
import random
from pathlib import Path
from translations import LANGUAGES
from core.chem_utils import safe_float, smiles_to_3d_block, get_pubchem_data, get_chembl_data, calculate_conformer_energies, compute_gasteiger_charges_block, get_quantum_descriptors, prepare_ligand_for_docking, run_ai_target_screening

# Временный try-except для импорта Тьютора
try:
    from core.bot import tutor_dialog
    # st.success("Бот импортирован")  
except Exception as e:
    st.error(f"Ошибка импорта бота: {e}")
        
st.set_page_config(page_title="BioSynth-EDU", layout="wide")

# --- 2. ИНИЦИАЛИЗАЦИЯ SESSION STATE ---
if 'mol_block' not in st.session_state:
    st.session_state.mol_block = None
if 'lang' not in st.session_state:
    st.session_state.lang = "Русский"

# --- 1. ФУНКЦИЯ ДЛЯ ОТОБРАЖЕНИЯ PDF ---
def display_pdf(url):
    # ссылка для предпросмотра
    embed_url = url.replace('/view?usp=drivesdk', '/preview')
    embed_url = embed_url.replace('/view', '/preview')
    
    # iframe для отображения
    pdf_display = f'<iframe src="{embed_url}" width="100%" height="600px" style="border:none;"></iframe>'
    st.markdown(pdf_display, unsafe_allow_html=True)

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

# --- 0. ФУНКЦИИ ВЗАИМНОГО СБРОСА (Колбэки, которые выполняются ДО рендера) ---
def on_kaz_change():
    # Если выбрали что-то в казахстанском каталоге — сбрасываем мировой и ручной ввод
    if st.session_state.kaz_select != t["select_placeholder"]:
        st.session_state.world_select = t["select_placeholder"]
        # Обновляем активный SMILES из казахстанского словаря
        st.session_state.active_smiles = kaz_options_global[st.session_state.kaz_select]

def on_world_change():
    # Если выбрали мировой пример — сбрасываем казахстанский каталог
    if st.session_state.world_select != t["select_placeholder"]:
        st.session_state.kaz_select = t["select_placeholder"]
        # Обновляем активный SMILES из мирового словаря
        st.session_state.active_smiles = world_examples_global[st.session_state.world_select]

# Инициализируем базовый SMILES в сессии, если приложения только открылось
if "active_smiles" not in st.session_state:
    st.session_state.active_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Аспирин по дефолту


# --- 4. БОКОВАЯ ПАНЕЛЬ ---
st.sidebar.header(t["sidebar_select_mol"])

# КАЗАХСТАНСКИЙ КАТАЛОГ
st.sidebar.subheader(t["sidebar_kaz_cat"])
kaz_options = {}
for m in catalog:
    display_name = m.get('name_local', {}).get(L_CODE, m.get('name', 'Unknown'))
    class_name = m.get('classification_local', {}).get(L_CODE, m.get('classification', 'Bioactiv'))
    kaz_options[f"{display_name} ({class_name})"] = m['smiles']

# Сохраняем в global для доступа внутри колбэка
kaz_options_global = kaz_options

selected_kaz = st.sidebar.selectbox(
    t["sidebar_kaz_label"], 
    options=[t["select_placeholder"]] + list(kaz_options.keys()),
    key="kaz_select",
    on_change=on_kaz_change  # Жесткий триггер сброса мирового каталога
)

# СТАНДАРТНЫЕ ПРИМЕРЫ
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
world_examples_global = world_examples

selected_world = st.sidebar.selectbox(
    t["sidebar_world_label"], 
    options=[t["select_placeholder"]] + list(world_examples.keys()),
    key="world_select",
    on_change=on_world_change  # Жесткий триггер сброса казахского каталога
)

# === УПРАВЛЕНИЕ РУЧНЫМ ВВОДОМ И ОЧИСТКОЙ ПАМЯТИ ===
st.sidebar.markdown("---")
st.sidebar.header(t["sidebar_manual"])

# Инициализируем базовый SMILES в сессии, если приложение только открылось
if "active_smiles" not in st.session_state:
    st.session_state.active_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Аспирин по дефолту

def on_manual_smiles_change():
    # Раз ввели руками или сработал триггер — сбрасываем оба селектбокса каталогов в плейсхолдер
    st.session_state.kaz_select = t["select_placeholder"]
    st.session_state.world_select = t["select_placeholder"]
    # Сбрасываем кэш старого докинга и 3D-модели
    st.session_state.prepared_pdbqt = None
    st.session_state.mol_block = None

# ЖЕСТКАЯ СВЯЗКА: делаем ключ текстового поля прямо нашей переменной сессии active_smiles.
# Убираем параметр value=..., чтобы Streamlit не затирал память дефолтами при перезапуске!
smiles = st.sidebar.text_input(
    t["sidebar_manual_label"], 
    key="active_smiles",
    on_change=on_manual_smiles_change
)

# === СИНХРОНИЗАЦИЯ ИНФОРМАЦИИ О ТЕКУЩЕЙ МОЛЕКУЛЕ (current_mol) ===
# Если выбран казахский каталог — вытаскиваем паспорт молекулы из базы данных
if st.session_state.kaz_select != t["select_placeholder"]:
    current_smiles = kaz_options[st.session_state.kaz_select]
    st.session_state.current_mol = next((m for m in catalog if m.get('smiles') == current_smiles), None)
else:
    st.session_state.current_mol = None

# Сохраняем локальное имя для бота/интерфейса
if st.session_state.current_mol:
    st.session_state.selected_mol_name = st.session_state.current_mol.get('name')
else:
    st.session_state.selected_mol_name = None

# 4. ОСНОВНОЙ ИНТЕРФЕЙС
st.title(f"🧪 {t.get('title_main', 'BioSynth-EDU')}")

tab_names = [
    t.get("tab_3d", "3D"), 
    t.get("tab_admet", "ADMET"), 
    t.get("tab_docking", "Docking"), 
    t.get("tab_edu", "Education"),
    t.get("tab_project", "Project") 
]
tab1, tab2, tab3, tab4, tab5 = st.tabs(tab_names)

with tab1:
    col1, col2 = st.columns([3, 1])
    with col1:
        c1, c2 = st.columns(2)
        if c1.button(t.get("btn_build_3d", "Построить 3D"), use_container_width=True):
            st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=False)
        
        if c2.button(t.get("btn_optimize", "Оптимизировать"), use_container_width=True):
            with st.spinner(t.get("spinner_optimize", "Расчет...")):
                st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=True)

        if st.session_state.mol_block:
            st.markdown(f"**{t.get('style_label', 'Стиль отображения')}:**")
            style_cols = st.columns(4)
            
            if 'viz_style' not in st.session_state:
                st.session_state.viz_style = 'stick'

            if style_cols[0].button("Stick", use_container_width=True): st.session_state.viz_style = 'stick'
            if style_cols[1].button("Sphere", use_container_width=True): st.session_state.viz_style = 'sphere'
            if style_cols[2].button("Line", use_container_width=True): st.session_state.viz_style = 'line'
            if style_cols[3].button("Surface", use_container_width=True): st.session_state.viz_style = 'surface'

            # --- АВТОМАТИЧЕСКИЙ БЛОК ОТРИСОВКИ 3D (ЗАВЯЗАН НА ИНДЕКС ВКЛАДКИ) ---
            view = py3Dmol.view(width=None, height=450)
            
            # Проверяем, открыта ли сейчас вкладка "Заряды" (ее индекс равен 1)
            is_charge_tab_active = st.session_state.get("srsp_tab_index") == 1
            
            if is_charge_tab_active and st.session_state.get('mol_block_charges'):
                # Рендерим молекулу с зарядами
                view.addModel(st.session_state.mol_block_charges, "mol")
                view.setColorByProperty({'prop': 'partialCharge', 'gradient': 'rwb', 'min': -0.3, 'max': 0.3})
            else:
                # В противном случае (вкладки 0, 2 или обычный режим) — стандартный вид
                if st.session_state.get('mol_block'):
                    view.addModel(st.session_state.mol_block, "mol")
                else:
                    active_sm = st.session_state.get('active_smiles', 'CC(=O)OC1=CC=CC=C1C(=O)O')
                    from core.chem_utils import smiles_to_3d_block
                    st.session_state.mol_block = smiles_to_3d_block(active_sm, optimize=True)
                    view.addModel(st.session_state.mol_block, "mol")
                
                # Применяем выбранный студентом стиль отображения структуры
                viz_st = st.session_state.get('viz_style', 'stick')
                if viz_st == 'stick':
                    view.setStyle({'stick': {'radius': 0.25}})
                elif viz_st == 'sphere':
                    view.setStyle({'sphere': {'scale': 0.9}})
                elif viz_st == 'line':
                    view.setStyle({'line': {'linewidth': 2}})
                elif viz_st == 'surface':
                    view.setStyle({'stick': {'radius': 0.1}})
                    view.addSurface(py3Dmol.VDW, {'opacity': 0.5, 'color': 'white'})
            
            view.zoomTo()
            view.setBackgroundColor('#ffffff')
            
            html_content = view._make_html().replace('width: 700px', 'width: 100%')
            components.html(html_content, height=480)
            
            try:
                prefix = selected_kaz.split()[0] if selected_kaz != t.get("select_placeholder") else "molecule"
            except:
                prefix = "molecule"
            
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            unique_filename = f"{prefix}_{timestamp}.sdf"

            st.download_button(
                label=t.get("download_sdf", "Скачать SDF"),
                data=st.session_state.mol_block,
                file_name=unique_filename,
                mime="chemical/x-mdl-sdfile",
                use_container_width=True
            )

   # --- ОБНОВЛЕННЫЙ БЛОК УПРАВЛЕНИЯ ЛАБОРАТОРНЫМИ (БЕЗ KX_MODE) ---
            st.divider()
            st.subheader("🔬 Квантово-химический экспресс-анализ (СРСП)")
            
            current_active_smiles = st.session_state.get("active_smiles", "")
            
            if not current_active_smiles:
                st.info("⚠️ Выберите молекулу для активации экспресс-анализа.")
            else:
                # Присваиваем вкладкам key, чтобы Streamlit сам отслеживал, где находится студент
                kx_tabs = st.tabs(["📊 Конформации", "⚡ Заряды (Эл. плотность)", "📝 Дескрипторы РС"], key="srsp_tab_index")
                
                # Вкладка 1: Конформации (Индекс 0)
                with kx_tabs[0]:
                    if st.button("Сгенерировать и рассчитать конформеры", use_container_width=True):
                        with st.spinner("Поиск глобального минимума..."):
                            energies, best_sdf = calculate_conformer_energies(current_active_smiles)
                            st.session_state.conf_energies = energies
                            st.session_state.best_sdf_block = best_sdf
                    
                    if 'conf_energies' in st.session_state and st.session_state.conf_energies:
                        chart_data = pd.DataFrame({
                            "Конформер №": [f"№{i+1}" for i in range(len(st.session_state.conf_energies))],
                            "Энергия (ккал/моль)": st.session_state.conf_energies
                        }).set_index("Конформер №")
                        st.bar_chart(chart_data)
                        
                        if st.session_state.get('best_sdf_block'):
                            st.download_button(
                                label="📥 Скачать самый стабильный конформер (SDF)",
                                data=st.session_state.best_sdf_block,
                                file_name="global_minimum_conformer.sdf",
                                mime="chemical/x-mdl-sdfile",
                                type="primary",
                                use_container_width=True
                            )
                
                # Вкладка 2: Заряды (Индекс 1)
                with kx_tabs[1]:
                    # Заряды считаются автоматически ОДИН раз при переходе на вкладку и кешируются в сессию
                    if 'last_charged_smiles' not in st.session_state or st.session_state.last_charged_smiles != current_active_smiles:
                        with st.spinner("Расчет парциальных зарядов по Гастейгеру..."):
                            res_charges = compute_gasteiger_charges_block(current_active_smiles)
                            if res_charges:
                                st.session_state.mol_block_charges = res_charges
                                st.session_state.last_charged_smiles = current_active_smiles
                                # Принудительно обновляем интерфейс, чтобы 3D-окно выше перекрасилось
                                st.rerun()
                    
                    st.success("✅ Электронная плотность успешно рассчитана и выведена на 3D-модель выше!")
                    st.caption("🔴 **Красный цвет:** дефицит электронов (электрофильный центр). 🔵 **Синий цвет:** избыток плотности (нуклеофильный центр).")
                
                # Вкладка 3: Дескрипторы (Индекс 2)
                with kx_tabs[2]:
                    kx_descriptors = get_quantum_descriptors(current_active_smiles)
                    if kx_descriptors is not None:
                        st.dataframe(kx_descriptors, use_container_width=True, hide_index=True)
                    else:
                        st.error("Ошибка дескрипторов. RDKit не смог верифицировать структуру SMILES.")
                        
        else:
            st.info(t.get("info_select_mol", "Выберите молекулу"))
            
    with col2:
        st.subheader(t.get("header_ref", "Свойства"))
        data = get_pubchem_data(smiles)
        
        if data:
            st.metric(t.get("metric_mw", "MW"), f"{data['mw']} g/mol")
            st.metric(t.get("metric_logp", "LogP"), data['logp'])
            st.metric(t.get("metric_rot_bonds", "Rot.Bonds"), data['rotatable_bonds'])
            
            st.divider()
            with st.spinner(t.get("spinner_chembl", "Поиск в ChEMBL...")):
                chembl_info = get_chembl_data(data['inchikey'])
            
            if chembl_info:
                st.write(f"🧬 **ChEMBL ID:** `{chembl_info['chembl_id']}`")
                phase = chembl_info['max_phase']
                status_color = "green" if phase == 4 else "orange"
                status_text = f"Phase {phase}"
                if phase == 4: status_text += f" ({t.get('status_approved', 'Approved')})"
                st.markdown(f"**{t.get('status_label', 'Статус:')}** <span style='color:{status_color}'>{status_text}</span>", unsafe_allow_html=True)
                
                with st.expander(t.get("mechanism_label", "Механизм действия")):
                    for m in chembl_info['mechanisms']:
                        st.write(f"• {m}")
            else:
                st.caption(t.get("no_chembl", "Данные ChEMBL не найдены"))

            st.divider()
            st.write(t.get("ext_links", "Внешние ссылки:"))
            
            inchikey = data['inchikey']
            st.link_button(t.get("btn_pubchem", "PubChem"), f"https://pubchem.ncbi.nlm.nih.gov/#query={inchikey}", use_container_width=True)
            st.link_button(t.get("btn_chembl", "ChEMBL"), f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={inchikey}", use_container_width=True)
            st.link_button(t.get("btn_similarity", "Similarity Search (70%)"), f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={smiles}&search_type=similarity&similarity=70", use_container_width=True, type="primary")
        else:
            st.warning(t.get("warn_no_pubchem", "Данные не найдены"))

#---Вкладка ADMET---
with tab2:
    st.header(t["admet_header"])
    st.markdown(t["admet_instructions"])
    
    col_btn1, col_btn2 = st.columns(2)
    with col_btn1:
        st.link_button(t["btn_open_admetlab"], "https://admetlab3.scbdd.com/", use_container_width=True)
    with col_btn2:
        st.link_button(t["btn_open_swissadme"], "http://www.swissadme.ch/", use_container_width=True)
    
    st.divider()
    
    uploaded_file = st.file_uploader(t["uploader_label"], type="csv")
    
    if uploaded_file:
        try:
            df = pd.read_csv(uploaded_file)
            
            with st.expander(t["expander_raw"]):
                st.dataframe(df)
                
            st.subheader(t["subheader_interp"])
            cols = df.columns.tolist()
            c1, c2, c3 = st.columns(3)

            logp_col = next((c for c in cols if 'logp' in c.lower()), None)
            if logp_col:
                val = safe_float(df[logp_col].iloc[0])
                with c1:
                    st.metric(t["metric_lipophilicity"], f"{val:.2f}")
                    if -1 < val < 5: 
                        st.success(t["status_optimal"])
                    else:
                        st.warning(t["status_extreme"])

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
            
            st.subheader(t.get("header_lipinski", "Соответствие правилам Drug-like"))
            
            violations = 0
            details = []

            mw_col = next((c for c in df.columns if 'mw' in c.lower() or 'weight' in c.lower()), None)
            logp_col = next((c for c in df.columns if 'logp' in c.lower()), None)
            hbd_col = next((c for c in df.columns if 'hbd' in c.lower() or 'donor' in c.lower()), None)
            hba_col = next((c for c in df.columns if 'hba' in c.lower() or 'acceptor' in c.lower()), None)

            if mw_col:
                val = safe_float(df[mw_col].iloc[0])
                if val > 500:
                    violations += 1
                    details.append(f"MW ({val:.1f}) > 500")

            if logp_col:
                val = safe_float(df[logp_col].iloc[0])
                if val > 5:
                    violations += 1
                    details.append(f"LogP ({val:.2f}) > 5")

            if hbd_col:
                val = safe_float(df[hbd_col].iloc[0])
                if val > 5:
                    violations += 1
                    details.append(f"HBD ({val}) > 5")

            if hba_col:
                val = safe_float(df[hba_col].iloc[0])
                if val > 10:
                    violations += 1
                    details.append(f"HBA ({val}) > 10")

            if violations == 0:
                st.balloons()
                st.success(t.get("lipinski_success", "✅ Молекула полностью соответствует правилу Липинского (Drug-like)."))
            else:
                st.warning(f"{t.get('lipinski_warn', 'Найдено нарушений: ')} {violations}.")
                for detail in details:
                    st.write(f"❌ {detail}")
                st.info(t.get("lipinski_info", "С точки зрения фармакокинетики, пероральный прием данного соединения может быть затруднен."))

        except Exception as e:
            st.error(f"{t.get('error_interp', 'Ошибка интерпретации: ')} {e}")
            
# --- Вкладка Докинг ---                   
with tab3:
    st.header("🧬 Молекулярный докинг: Подбор сайта связывания")
    
    # ЭТАП 1: Подготовка лиганда
    st.header(t["docking_header"])
    
    if st.session_state.get('mol_block') or smiles:
        st.success(t["docking_mol_ready"])
        
        col_prep1, col_prep2 = st.columns(2)
        
        with col_prep1:
            st.markdown(t["docking_checklist_title"])
            st.markdown(t["docking_checklist_items"])
            
            if st.button(t["btn_run_prep"], use_container_width=True):
                with st.spinner(t["spinner_meeko"]):
                    pdbqt_data = prepare_ligand_for_docking(smiles)
                    if pdbqt_data:
                        st.session_state.prepared_pdbqt = pdbqt_data
                        st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=True)
                        st.balloons()
                    else:
                        st.error(t["docking_error"])

        with col_prep2:
            st.info(t["docking_note_student"])
            
            if st.session_state.get('prepared_pdbqt'):
                st.download_button(
                    label=t["btn_download_pdbqt"],
                    data=st.session_state.prepared_pdbqt,
                    file_name="ligand.pdbqt",
                    mime="text/plain",
                    use_container_width=True
                )
                
                st.write(t["docking_view_label"])
                
                try:
                    view = py3Dmol.view(width=400, height=400)
                    view.addModel(st.session_state.mol_block, "sdf")
                    view.setStyle({'stick': {'color': 'spectrum', 'radius': 0.15}, 'sphere': {'scale': 0.25}})
                    view.zoomTo()
                    view.setBackgroundColor('#ffffff')
                    components.html(view._make_html(), height=410)
                    st.caption(t["docking_caption"])
                except Exception as e:
                    st.error(f"Render error: {e}")
            
        st.divider()
        st.subheader(t["docking_next_steps_header"])
        st.write(t["docking_next_steps_text"])

        # Векторный скрининг и выбор биологической мишени scPDB
        st.subheader(t["docking_stage2_title"])
        st.write(t["docking_stage2_desc"])
        
        if st.button(t["btn_run_screening"], use_container_width=True):
            with st.spinner(t["spinner_screening"]):
                res = run_ai_target_screening(smiles)
                
            if "error" in res:
                st.error(f"Error: {res['error']}")
            else:
                st.success(t["screening_success"])
                
                top = res["top_match"]
                
                st.info(
                    t["top_match_title"].format(pdb_id=top['id']) + "\n\n" +
                    t["top_match_reason"].format(sim=top['sim']) + "\n\n" +
                    t["top_match_instruction"].format(pdb_id=top['id'])
                )

                df_data = []
                for item in res["all_candidates"]:
                    df_data.append({
                        t["col_pdb_id"]: item["id"],
                        t["col_tanimoto"]: f"{item['sim']:.4f}",
                        t["col_pkd"]: f"{item['score']:.2f}"
                    })
                
                df = pd.DataFrame(df_data)
                
                st.write(t["table_title"])
                st.dataframe(df, use_container_width=True, hide_index=True)
                               
                if top['id']:
                    # Вырезаем первые 4 символа (PDB ID) из ID сайта scPDB
                    best_pdb = top['id'][:4].upper()
                    pdb_url = f"https://www.rcsb.org/structure/{best_pdb}"
                    
                    st.write("") # Небольшой отступ
                    st.link_button(
                        label=t["btn_go_to_pdb"].format(pdb_id=best_pdb),
                        url=pdb_url,
                        use_container_width=True,
                        type="primary" # Сделаем её яркой и акцентной
                    )
        
# --- Вкладка Обучение ---
with tab4:
    current_mol = next((m for m in catalog if m['smiles'] == smiles), None)
    if current_mol:
        with st.expander("🇰🇿 Сведения о казахстанской разработке", expanded=True):
            st.markdown(f"### Препарат: {current_mol['name']}")
            # Используем данные из вашего профиля о казахстанских ученых
            st.write(f"**Авторы:** {', '.join(current_mol.get('authors', []))}")
            st.info(f"**Краткое описание:** {current_mol.get('description', '')}")
        st.markdown("---") 

    itabs = st.tabs(["📖 Пособие", "🎥 Видео", "🧪 Лабораторные", "🔬 Магистрантам"])

    # РАЗДЕЛ 1: ПОСОБИЕ 
    with itabs[0]:
        st.subheader("Учебные пособия")
        
        # Выбор пособия
        manual_choice = st.radio(
            "Выберите литературу для изучения:",
            ["Компьютерная химия: прогнозирование", "Квантово-химические расчеты и молекулярное моделирование"],
            horizontal=True
        )

        # Словарь ссылок
        manual_links = {
            "Компьютерная химия: прогнозирование": "https://drive.google.com/file/d/1CxTdwJmiaqsa54wZoODhADYIqHoIk23g/view?usp=drivesdk",
            "Квантово-химические расчеты и молекулярное моделирование": "https://drive.google.com/file/d/1JeO18y1QNPl-YzA_aRTSe8owaNMvms7k/view?usp=drivesdk"
        }

        col_pdf, col_tutor = st.columns([2, 1])

        with col_pdf:
            # Отображаем выбранный PDF
            display_pdf(manual_links[manual_choice])
            st.caption(f"Просмотр файла: {manual_choice}. Если файл не отображается, проверьте доступ по ссылке.")

        with col_tutor:
            st.subheader("🤖 ИИ-Тьютор")
            st.write("Используйте материалы пособия для подготовки вопросов.")
            
            #  ИИ-Тьютор 
            if st.button("💬 Открыть диалог с Тьютором", type="primary", use_container_width=True):
                tutor_dialog()

    # Остальные разделы
    with itabs[1]:
        st.video("https://www.youtube.com/watch?v=dQw4w9WgXcQ")
        st.radio("Тест: Как изменится LogP?", ["Увеличится", "Уменьшится"], key="v_test")

    with itabs[2]:
        st.info("Лабораторные работы по цифровой трансформации химии") #
    with itabs[3]:
        st.info("Кейсы")
        
# --- Вкладка 5: ИССЛЕДОВАТЕЛЬСКИЙ ПРОЕКТ ---
with tab5:
    st.header(t.get("tab_project", "🚀 Исследовательский проект"))

    st.markdown(t.get(
        "project_desc",
        "**Исследовательский проект по органической химии для студентов и учащихся**"
    ))

    st.markdown(t.get(
        "project_start_hint",
        "👈 Для начала работы выберите молекулу в боковой панели"
    ))

    st.divider()

    mol_data = st.session_state.get('current_mol')
    project_smiles = mol_data.get('smiles', smiles) if mol_data else smiles

    if project_smiles and project_smiles.strip():

        if mol_data:
            st.subheader(f"🔬 {mol_data.get('name', t.get('selected_mol', 'Выбранное соединение'))}")
            if mol_data.get('name_local'):
                st.caption(mol_data.get('name_local', {}).get(L_CODE, ''))
        else:
            st.subheader(t.get("input_struct", "🔬 Введённая структура"))

        st.info(f"**{t.get('current_smiles', 'Текущий SMILES')}:** `{project_smiles}`")

        if mol_data:
            st.markdown(f"### {t.get('kz_info', '🇰🇿 Сведения о разработке')}")
            col_info1, col_info2 = st.columns([2, 1])

            with col_info1:
                st.write(f"**{t.get('authors', 'Авторы')}:** {', '.join(mol_data.get('authors', ['—']))}")
                st.info(f"**{t.get('description', 'Описание')}:** {mol_data.get('description', t.get('info_missing', '—'))}")

            with col_info2:
                if mol_data.get('classification'):
                    st.metric(t.get('classification', 'Классификация'), mol_data.get('classification', '—'))
                if mol_data.get('year'):
                    st.metric(t.get('year', 'Год'), mol_data.get('year', '—'))
            st.divider()

        st.subheader(t.get("task_title", "📝 Задание на проект"))
        with st.expander(t.get("task_full", "Открыть полное задание"), expanded=True):
            st.markdown(t.get("task_step_1", "1. Прогноз PASS Online."))
            st.markdown(t.get("task_step_2", "2. Оценка ADMET и правила Липинского."))
            st.markdown(t.get("task_step_3", "3. Модификация структуры в редакторе."))
            st.markdown(t.get("task_step_4", "4. Выводы для доклада."))
        st.divider()

        st.subheader(t.get("editor_title", "🧪 Редактор структуры"))
        st.markdown(f"**{t.get('editor_task', 'Задание: Измените структуру и примените изменения.')}**")

        editor_key = f"ketcher_{hash(project_smiles)}_{st.session_state.lang}"
        edited = st_ketcher(project_smiles, key=editor_key)

        if edited and edited != project_smiles:
            st.success(t.get("structure_changed", "✅ Структура изменена!"))
            col_btn1, col_btn2 = st.columns(2)

            with col_btn1:
                if st.button(t.get("apply_btn", "🔄 Применить изменения"), use_container_width=True, type="primary"):
                    if mol_data:
                        st.session_state.current_mol['smiles'] = edited
                    st.session_state.active_smiles = edited
                    st.rerun()

            with col_btn2:
                st.download_button(
                    label=t.get("download_btn", "💾 Скачать SMILES"),
                    data=edited,
                    file_name=f"modified_{datetime.datetime.now().strftime('%H%M')}.smi",
                    mime="text/plain",
                    use_container_width=True
                )
        else:
            st.caption(t.get("change_editor_hint", "👆 Измените молекулу в редакторе выше"))

        st.divider()

        st.subheader("📚 " + t.get("tools_header", "Обучение"))
        
        l_col1, l_col2, l_col3 = st.columns(3)
        l_col1.link_button(t.get("pass_lecture", "🎥 PASS"), "https://drive.google.com/file/d/1nqHlBMpj6RZ28MvSElHiHgt495GCcSwT/view?usp=drive_link", use_container_width=True)
        l_col2.link_button(t.get("adme_lecture", "🎥 ADME"), "https://drive.google.com/file/d/1WggF966FXCrgO41QfJxtU9Ls8MgZjJbB/view?usp=drive_link", use_container_width=True)
        l_col3.link_button(t.get("pubchem_lecture", "🎥 PubChem"), "https://drive.google.com/file/d/1MTygaMnjEcuIzL0cSOV391QkwBUQcUO7/view?usp=sharing", use_container_width=True)

        s_col1, s_col2, s_col3 = st.columns(3)
        s_col1.link_button("🌐 PASS Online", "http://www.way2drug.com/passonline/", use_container_width=True, type="primary")
        s_col2.link_button("🧪 SwissADME", "http://www.swissadme.ch/", use_container_width=True, type="primary")
        s_col3.link_button("📊 PubChem", f"https://pubchem.ncbi.nlm.nih.gov/#query={project_smiles}", use_container_width=True, type="primary")

        st.divider()
        if st.button(t.get("quiz_btn", "📝 Пройти тест"), use_container_width=True, type="primary"):
            from core.assessment import get_assessment_data, run_quiz_dialog
            
            tests, open_qs, cols_map = get_assessment_data()
            if tests:
                run_quiz_dialog(tests, open_qs, cols_map, t)
            else:
                st.error("Не удалось загрузить данные тестов.")

    else:
        st.warning(t.get("mol_not_selected", "⚠️ Молекула не выбрана"))
        st.info(f"**{t.get('how_start', 'Как начать:')}**\n{t.get('start_step_1', '• Выберите соединение или введите SMILES в боковой панели')}")

# Кнопка Тьютора в сайдбаре
st.sidebar.divider()
if st.sidebar.button("💬 Задать вопрос Тьютору", use_container_width=True, type="primary"):
    tutor_dialog()
