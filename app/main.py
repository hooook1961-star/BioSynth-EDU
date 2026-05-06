import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import pandas as pd
from core.chem_utils import smiles_to_3d_block, get_pubchem_data

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

tab1, tab2, tab3 = st.tabs(["🔬 3D Структура и Кванты", "📊 ADMET Анализ", "📖 Обучение"])

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
        st.subheader("Характеристики")
        data = get_pubchem_data(smiles)
        if data:
            st.metric("М. вес", f"{data['mw']} г/моль")
            st.metric("LogP", data['logp'])
            st.write(f"**Формула:** {data['formula']}")
            st.divider()
            st.caption("Данные подгружены из PubChem")

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
    st.header("Курс лекций и тесты")
    st.video("https://www.youtube.com/watch?v=dQw4w9WgXcQ")
    st.markdown("### Мини-тест\nКак изменится липофильность (LogP) при добавлении -OH группы?")
    st.radio("Выберите ответ:", ["Увеличится", "Уменьшится", "Не изменится"])
