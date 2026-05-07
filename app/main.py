import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import pandas as pd
import json
import os
import datetime
from core.chem_utils import smiles_to_3d_block, get_pubchem_data, get_chembl_data, prepare_ligand_for_docking

# --- БЛОК ЗАГРУЗКИ ДАННЫХ ---
# Выходим из папки app в корень и заходим в data
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(BASE_DIR, 'data', 'catalog.json')

@st.cache_data
def load_catalog():
    try:
        with open(DATA_PATH, 'r', encoding='utf-8') as f:
            return json.load(f)
    except FileNotFoundError:
        # Если файла нет, возвращаем пустой список, чтобы приложение не упало
        return []

catalog = load_catalog()
# Создаем словарь для боковой панели: { "Название (Класс)": "SMILES" }
examples = {f"{m['name']} ({m.get('classification', 'Биоактив')})": m['smiles'] for m in catalog}
# ----------------------------------

# 1. КОНФИГУРАЦИЯ
st.set_page_config(page_title="BioSynth-EDU", layout="wide")

# 2. ИНИЦИАЛИЗАЦИЯ СОСТОЯНИЯ
if 'mol_block' not in st.session_state:
    st.session_state.mol_block = None

# 3. БОКОВАЯ ПАНЕЛЬ
st.sidebar.header("🧪 Выбор молекулы")

# --- ГРУППА 1: КАЗАХСТАНСКИЙ КАТАЛОГ (BioSynth-EDU) ---
st.sidebar.subheader("🇰🇿 Разработки Казахстана")
# Сохраняем названия и действия из JSON
kaz_options = {f"{m['name']} ({m.get('classification', 'Биоактив')})": m['smiles'] for m in catalog}

selected_kaz = st.sidebar.selectbox(
    "Отечественные препараты и кейсы:", 
    options=["-- Выберите из списка --"] + list(kaz_options.keys())
)

# --- ГРУППА 2: ПРИМЕРЫ ЛЕКАРСТВ ---
st.sidebar.subheader("🌍 Стандартные примеры")
examples = {
    "Аспирин (Анальгетик)": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Кофеин (Стимулятор)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Парацетамол (Жаропонижающее)": "CC(=O)NC1=CC=C(O)C=C1",
    "Ибупрофен (НПВС)": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Пенициллин G (Антибиотик)": "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
    "Никотин (Алкалоид)": "CN1CCCC1C2=CN=CC=C2",
    "Дофамин (Нейромедиатор)": "C1=CC(=C(C=C1CCN)O)O"
}

selected_world = st.sidebar.selectbox(
    "Примеры лекарственных веществ:", 
    options=["-- Выберите из списка --"] + list(examples.keys())
)

st.sidebar.markdown("---")

# --- ЛОГИКА ОПРЕДЕЛЕНИЯ ТЕКУЩЕГО SMILES ---
# По умолчанию ставим первый элемент из мировых примеров (Аспирин)
current_smiles = examples["Аспирин (Анальгетик)"]

# Если выбран Казахстанский препарат - берем его
if selected_kaz != "-- Выберите из списка --":
    current_smiles = kaz_options[selected_kaz]
# Если выбран Мировой пример - берем его
elif selected_world != "-- Выберите из списка --":
    current_smiles = examples[selected_world]

# --- ПОЛЕ ВВОДА SMILES (ДЛЯ РЕДАКТИРОВАНИЯ) ---
st.sidebar.header("✍️ Ввести SMILES")
# Поле text_input динамически получает smiles из выбранного списка
smiles = st.sidebar.text_input("Или вставьте SMILES ниже:", value=current_smiles)

# 4. ОСНОВНОЙ ИНТЕРФЕЙС
st.title("🧪 BioSynth-EDU: Исследовательская платформа")

tab1, tab2, tab3, tab4 = st.tabs(["🔬 3D Структура", "📊 ADMET Анализ", "🧬 Докинг", "📖 Обучение"])

with tab1:
    col1, col2 = st.columns([3, 1])
    
    with col1:
        c1, c2 = st.columns(2)
        if c1.button("🏗️ Построить 3D", use_container_width=True):
            st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=False)
        
        if c2.button("✨ Оптимизировать (MMFF94)", use_container_width=True):
            with st.spinner("Минимизация энергии..."):
                st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=True)

        if st.session_state.mol_block:
            # Визуализация
            view = py3Dmol.view(width=700, height=500)
            view.addModel(st.session_state.mol_block, "mol")
            view.setStyle({'stick': {'radius':0.2}, 'sphere': {'scale':0.3}})
            view.zoomTo()
            view.setBackgroundColor('#ffffff')
            components.html(view._make_html(), height=550)
            # 1. Формируем безопасное и уникальное имя файла
            # Пытаемся взять имя из выбора, если нет — используем "molecule"
            try:
                # Берем первое слово из названия (например, "Просидол")
                prefix = selected_name.split()[0]
            except (NameError, IndexError, AttributeError):
                prefix = "molecule"
            
            # Добавляем дату и время: ГГГГММДД_ЧЧММСС
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            unique_filename = f"{prefix}_{timestamp}.sdf"

            # 2. Кнопка скачивания SDF
            st.download_button(
                label="💾 Скачать структуру (SDF)",
                data=st.session_state.mol_block,
                file_name=unique_filename,
                mime="chemical/x-mdl-sdfile",
                help="SDF файл содержит 3D координаты атомов для работы в проф. софте (AutoDock, PyMOL и др.)"
            )
        else:
            st.info("Выберите молекулу слева и нажмите 'Построить 3D'")

with col2:
        st.subheader("📚 Справочник")
        data = get_pubchem_data(smiles)
        
        if data:
            st.metric("М. вес", f"{data['mw']} г/моль")
            st.metric("LogP", data['logp'])
            st.metric("Вращающихся связей", data['rotatable_bonds'])
            
            # --- НОВЫЙ БЛОК: Данные ChEMBL ---
            st.divider()
            with st.spinner("Запрос к ChEMBL..."):
                chembl_info = get_chembl_data(data['inchikey'])
            
            if chembl_info:
                st.write(f"🧬 **ChEMBL ID:** `{chembl_info['chembl_id']}`")
                
                # Статус одобрения
                phase = chembl_info['max_phase']
                status_color = "green" if phase == 4 else "orange"
                st.markdown(f"**Статус:** <span style='color:{status_color}'>Phase {phase} (Одобрено)</span>" if phase == 4 else f"**Статус:** Phase {phase}", unsafe_allow_html=True)
                
                # Механизмы
                with st.expander("🔬 Механизм действия"):
                    for m in chembl_info['mechanisms']:
                        st.write(f"• {m}")
            else:
                st.caption("Биологическая активность в ChEMBL не найдена")

            st.divider()
            
            # Блок внешних ссылок 
            st.write("🔗 **Внешние базы:**")
            
            # Ссылка на PubChem
            pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/#query={data['inchikey']}"
            st.link_button("Профиль в PubChem", pubchem_url, use_container_width=True)
            
            # Ссылка на ChEMBL 
            chembl_url = f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={data['inchikey']}"
            st.link_button("Данные ChEMBL (IC50/Ki)", chembl_url, use_container_width=True)
            
            # Кнопка сходства
            chembl_sim_url = f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={smiles}&search_type=similarity&similarity=70"
            st.link_button("Найти похожие", chembl_sim_url, use_container_width=True, type="primary")

        else:
            st.warning("Данные в PubChem не найдены")

# Вкладка ADMET 
with tab2:
    st.header("📊 ADMET Анализ и интерпретация")
    
    # Блок с инструкциями и кнопками
    st.markdown("""
    Для проведения анализа ADMET:
    1. Нажмите кнопку **«Открыть ADMETlab 3.0»** ниже.
    2. Скопируйте ваш SMILES и вставьте его в поле ввода на сайте.
    3. После завершения расчета скачайте результат в формате **CSV**.
    4. Загрузите файл с помощью кнопки Upload для автоматической интерпретации.
    """)
    
    # Кнопки со ссылками
    col_btn1, col_btn2 = st.columns(2)
    with col_btn1:
        st.link_button("🌐 Открыть ADMETlab 3.0", "https://admetlab3.scbdd.com/", use_container_width=True)
    with col_btn2:
        st.link_button("🧪 Альтернатива: SwissADME", "http://www.swissadme.ch/", use_container_width=True)
    
    st.divider()
    
    # Загрузка файла
    uploaded_file = st.file_uploader("📥 Загрузите CSV файл с результатами ADMETlab", type="csv")
    
    if uploaded_file:
        try:
            df = pd.read_csv(uploaded_file)
            
            with st.expander("👁️ Посмотреть все сырые данные CSV"):
                st.dataframe(df)
                
            st.subheader("📝 Краткая интерпретация ключевых параметров")
            cols = df.columns.tolist()
            c1, c2, c3 = st.columns(3)

            # Вспомогательная функция для безопасного превращения в число
            def safe_float(val):
                if isinstance(val, (int, float)): return float(val)
                val_str = str(val).strip().lower()
                if val_str == 'yes': return 1.0
                if val_str == 'no': return 0.0
                try: return float(val_str)
                except: return 0.0

            # 1. Липофильность (LogP)
            logp_col = next((c for c in cols if 'logp' in c.lower()), None)
            if logp_col:
                val = safe_float(df[logp_col].iloc[0])
                with c1:
                    st.metric("LogP (Липофильность)", f"{val:.2f}")
                    if -1 < val < 5: 
                        st.success("✅ Оптимально.")
                    else:
                        st.warning("⚠️ Крайние значения.")

            # 2. Гематоэнцефалический барьер (BBB)
            bbb_col = next((c for c in cols if 'bbb' in c.lower()), None)
            if bbb_col:
                raw_val = df[bbb_col].iloc[0]
                val = safe_float(raw_val)
                with c2:
                    st.metric("BBB (Проницаемость)", str(raw_val))
                    if val > 0.5 or str(raw_val).lower() == 'yes':
                        st.warning("🧠 Проникает через ГЭБ.")
                    else:
                        st.success("🛡️ Низкий риск для ЦНС.")

            # 3. Токсичность (hERG / Toxicity)
            herg_col = next((c for c in cols if 'herg' in c.lower() or 'tox' in c.lower()), None)
            if herg_col:
                raw_val = df[herg_col].iloc[0]
                val = safe_float(raw_val)
                with c3:
                    st.metric("Токсичность / hERG", str(raw_val))
                    if val > 0.5 or str(raw_val).lower() == 'yes':
                        st.error("💔 Высокий риск.")
                    else:
                        st.success("❤️ Риск низкий.")

            st.divider()
            
            # --- Анализ правил Липинского ---
            st.subheader("🧐 Соответствие правилам Drug-like")
            
            # Ищем колонки веса, доноров и акцепторов
            mw_col = next((c for c in cols if any(w in c.lower() for w in ['mw', 'weight'])), None)
            hbd_col = next((c for c in cols if 'hbd' in c.lower() or 'donor' in c.lower()), None)
            hba_col = next((c for c in cols if 'hba' in c.lower() or 'acceptor' in c.lower()), None)
            
            violations = 0
            if mw_col and safe_float(df[mw_col].iloc[0]) > 500: violations += 1
            if logp_col and safe_float(df[logp_col].iloc[0]) > 5: violations += 1
            if hbd_col and safe_float(df[hbd_col].iloc[0]) > 5: violations += 1
            if hba_col and safe_float(df[hba_col].iloc[0]) > 10: violations += 1
            
            if violations == 0:
                st.balloons()
                st.success("🌟 Молекула полностью соответствует правилу Липинского!")
            else:
                st.warning(f"⚠️ Нарушений правила Липинского: {violations}.")
                st.info("💡 Напоминание: Допускается 1 нарушение для сохранения 'drug-likeness'.")

        except Exception as e:
            st.error(f"❌ Ошибка интерпретации: {e}")
with tab3:
    st.header("🛠️ Подготовка лиганда к докингу")
    
    if st.session_state.mol_block:
        st.success("✅ 3D-структура обнаружена и готова к обработке.")
        
        col_prep1, col_prep2 = st.columns(2)
        
        with col_prep1:
            st.markdown("""
            **Чек-лист подготовки:**
            1. Добавление неявных водородов.
            2. Генерация 3D-конформации.
            3. Минимизация энергии (силовое поле MMFF94).
            4. **Определение активных торсионных углов (PDBQT).**
            """)
            
            if st.button("⚙️ Запустить полную подготовку", use_container_width=True):
                with st.spinner("Работают Meeko и RDKit: расчет зарядов и торсионов..."):
                    # Вызываем профессиональную подготовку
                    pdbqt_data = prepare_ligand_for_docking(smiles)
                    if pdbqt_data:
                        st.session_state.prepared_pdbqt = pdbqt_data
                        st.balloons()
                        st.info("Лиганд готов! Рассчитаны торсионы и заряды.")
                    else:
                        st.error("Ошибка при подготовке PDBQT.")

        with col_prep2:
            st.info("ℹ️ **Заметка для студентов:** Докинг имитирует 'ключ и замок'. Чтобы ключ подошел, он должен иметь правильные углы связей.")
            
            if 'prepared_pdbqt' in st.session_state:
                st.download_button(
                    label="📥 Скачать готовый PDBQT",
                    data=st.session_state.prepared_pdbqt,
                    file_name="ligand.pdbqt",
                    mime="text/plain",
                    use_container_width=True
                )
                
# --- БЛОК Визуализации подготовленного лиганда ---
                st.write("🔍 **Просмотр подготовленной структуры:**")
                import py3Dmol
                
                # Используем SDF блок, так как он лучше отображает связи в py3Dmol
                # Если SDF хранится в mol_block, берем его. Если нет - оставляем pdbqt, но с фиксом.
                mol_data = st.session_state.get('mol_block', st.session_state.prepared_pdbqt)
                mol_format = 'sdf' if st.session_state.get('mol_block') else 'pdbqt'

                view = py3Dmol.view(width=400, height=400) 
                view.addModel(mol_data, mol_format)
                
                # Добавляем настройку для корректного отображения связей
                view.setStyle({'stick': {'color': 'spectrum', 'radius': 0.15}, 'sphere': {'scale': 0.25}})
                view.zoomTo()
                
                # Рендерим компонент в Streamlit
                # Используем полное название, так как import streamlit.components.v1 as components в начале файла
                import streamlit.components.v1 as components
                components.html(view._make_html(), height=410)
                
                st.caption("Оптимизированная 3D-модель (подготовлено для анализа).")
            
        st.divider()
        st.subheader("🎓 Что дальше?")
        st.write("""
        После подготовки лиганда вам необходимо:
        1. Подготовить **белок-мишень** (удалить воду, добавить заряды).
        2. Определить **Grid Box** (координаты активного центра).
        3. Запустить расчет в AutoDock Vina или аналогичном ПО.
        """)
        
    else:
        st.warning("⚠️ Сначала постройте 3D модель на первой вкладке!")
        
with tab4:
    # --- БЛОК JSON ---
    current_mol = next((m for m in catalog if m['smiles'] == smiles), None)
    
    if current_mol:
        # данные из каталога
        with st.expander("🇰🇿 Сведения о казахстанской разработке", expanded=True):
            st.markdown(f"### Препарат: {current_mol['name']}")
            st.write(f"**Авторы:** {', '.join(current_mol.get('authors', []))}")
            st.write(f"**Патент:** {current_mol.get('patent', '—')}")
            st.info(f"**Краткое описание:** {current_mol.get('description', '')}")
        st.markdown("---") # Разделительная линия перед основным текстом
    # -------------------------------------------------------------

    st.write("Инструкции для студентов...")
    st.header("Курс лекций и тесты")
    st.video("https://www.youtube.com/watch?v=dQw4w9WgXcQ")
    st.markdown("### Мини-тест\nКак изменится липофильность (LogP) при добавлении -OH группы?")
    st.radio("Выберите ответ:", ["Увеличится", "Уменьшится", "Не изменится"])
