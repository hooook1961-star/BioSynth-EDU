import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import pandas as pd
from core.chem_utils import smiles_to_3d_block, get_pubchem_data, get_chembl_data

# 1. КОНФИГУРАЦИЯ
st.set_page_config(page_title="BioSynth-EDU", layout="wide")

# 2. ИНИЦИАЛИЗАЦИЯ СОСТОЯНИЯ
if 'mol_block' not in st.session_state:
    st.session_state.mol_block = None

# 3. БОКОВАЯ ПАНЕЛЬ
st.sidebar.header("🧪 Выберите молекулу")
examples = {
    "Аспирин (Анальгетик)": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Кофеин (Стимулятор)": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Парацетамол (Жаропонижающее)": "CC(=O)NC1=CC=C(O)C=C1",
    "Ибупрофен (НПВС)": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Пенициллин G (Антибиотик)": "CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
    "Никотин (Алкалоид)": "CN1CCCC1C2=CN=CC=C2",
    "Дофамин (Нейромедиатор)": "C1=CC(=C(C=C1CCN)O)O"
}

selected_name = st.sidebar.selectbox("Готовые примеры лекарственных веществ:", list(examples.keys()))

st.sidebar.markdown("---")
st.sidebar.header("✍️ Свой ввод")
smiles = st.sidebar.text_input("Или вставьте SMILES ниже:", examples[selected_name])

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
            
            # Кнопка скачивания SDF
            st.download_button(
                label="💾 Скачать структуру (SDF)",
                data=st.session_state.mol_block,
                file_name=f"{selected_name.split()[0]}.sdf",
                mime="chemical/x-mdl-sdfile",
                help="SDF файл содержит 3D координаты атомов для работы в проф. софте"
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
            
            # Блок внешних ссылок (теперь более точный)
            st.write("🔗 **Внешние базы:**")
            
            # Ссылка на PubChem
            pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/#query={data['inchikey']}"
            st.link_button("Профиль в PubChem", pubchem_url, use_container_width=True)
            
            # Ссылка на ChEMBL (теперь используем InChIKey для точности)
            chembl_url = f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={data['inchikey']}"
            st.link_button("Данные ChEMBL (IC50/Ki)", chembl_url, use_container_width=True)
            
            # Кнопка сходства
            chembl_sim_url = f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={smiles}&search_type=similarity&similarity=70"
            st.link_button("Найти похожие", chembl_sim_url, use_container_width=True, type="primary")

        else:
            st.warning("Данные в PubChem не найдены")

# Вкладки ADMET и Обучение остаются такими же, как мы обсуждали ранее
with tab2:
    st.header("Прогнозирование свойств (SwissADME / ADMETlab)")
    st.info("Используйте внешние сервисы для глубокого анализа, затем загрузите CSV результат.")
    
    serv_col1, serv_col2 = st.columns(2)
    with serv_col1:
        st.link_button("Открыть SwissADME ↗️", "http://www.swissadme.ch/", use_container_width=True)
    with serv_col2:
        st.link_button("Открыть ADMETlab 3.0 ↗️", "https://admetmesh.scbdd.com/", use_container_width=True)

    st.divider()
    uploaded_file = st.file_uploader("Загрузите CSV результат", type="csv")
    if uploaded_file:
        try:
            df = pd.read_csv(uploaded_file, sep=None, engine='python')
            st.success("Данные загружены!")
            st.dataframe(df.T, use_container_width=True)
        except Exception as e:
            st.error(f"Ошибка: {e}")
with tab3:
    st.header("🛠️ Подготовка лиганда к докингу")
    
    if st.session_state.mol_block:
        st.success("✅ 3D-структура обнаружена и готова к обработке.")
        
        col_prep1, col_prep2 = st.columns(2)
        
        with col_prep1:
            st.markdown("""
            **Чек-лист подготовки:**
            1.  Добавление неявных водородов (H-atoms).
            2.  Генерация 3D-конформации.
            3.  Минимизация энергии (силовое поле MMFF94).
            4.  Определение активных торсионных углов.
            """)
            
            if st.button("⚙️ Запустить полную подготовку", use_container_width=True):
                with st.spinner("Оптимизация геометрии..."):
                    # Принудительно оптимизируем для докинга
                    st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=True)
                st.balloons()
                st.info("Лиганд оптимизирован. Теперь он находится в локальном минимуме энергии.")

        with col_prep2:
            st.info("ℹ️ **Заметка для студентов:** Докинг имитирует 'ключ и замок'. Чтобы ключ подошел, он должен иметь правильные углы связей.")
            
          if 'prepared_pdbqt' in st.session_state:
            st.download_button(
                label="📥 Скачать готовый PDBQT",
                data=st.session_state.prepared_pdbqt,
                file_name="ligand.pdbqt",
                mime="text/plain"
            
        st.divider()
        st.subheader("🎓 Что дальше?")
        st.write("""
        После подготовки лиганда вам необходимо:
        1.  Подготовить **белок-мишень** (удалить воду, добавить заряды).
        2.  Определить **Grid Box** (координаты активного центра).
        3.  Запустить расчет в AutoDock Vina или аналогичном ПО.
        """)
        
    else:
        st.warning("⚠️ Сначала постройте 3D модель на первой вкладке!")
        
with tab4:
    st.header("Курс лекций и тесты")
    st.video("https://www.youtube.com/watch?v=dQw4w9WgXcQ")
    st.markdown("### Мини-тест\nКак изменится липофильность (LogP) при добавлении -OH группы?")
    st.radio("Выберите ответ:", ["Увеличится", "Уменьшится", "Не изменится"])
