import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
from core.chem_utils import smiles_to_3d_block, get_pubchem_data

# 1. КОНФИГУРАЦИЯ (Всегда первая!)
st.set_page_config(page_title="BioSynth-EDU", layout="wide")

# 2. ИНИЦИАЛИЗАЦИЯ СОСТОЯНИЯ (Чтобы данные не пропадали при нажатии кнопок)
if 'mol_block' not in st.session_state:
    st.session_state.mol_block = None

# 3. БОКОВАЯ ПАНЕЛЬ
st.sidebar.header("⚙️ Параметры")
examples = {
    "Фенол": "c1ccccc1O",
    "Аспирин": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Кофеин": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}
selected_name = st.sidebar.selectbox("Примеры", list(examples.keys()))
smiles = st.sidebar.text_input("SMILES", examples[selected_name])

# 4. ОСНОВНОЙ ИНТЕРФЕЙС
st.title("🧪 BioSynth-EDU")

tab1, tab2, tab3 = st.tabs(["🔬 Визуализация и Кванты", "📊 ADMET Анализ", "📖 Обучение"])

with tab1:
    col1, col2 = st.columns([3, 1])
    
    with col1:
        # Кнопки управления
        c1, c2 = st.columns(2)
        if c1.button("🏗️ Построить 3D"):
            st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=False)
        
        if c2.button("✨ Оптимизировать (MMFF94)"):
            with st.spinner("Минимизация энергии..."):
                st.session_state.mol_block = smiles_to_3d_block(smiles, optimize=True)

        # Визуализация
        if st.session_state.mol_block:
            view = py3Dmol.view(width=700, height=500)
            view.addModel(st.session_state.mol_block, "mol")
            view.setStyle({'stick': {'radius':0.2}, 'sphere': {'scale':0.3}})
            view.zoomTo()
            components.html(view._make_html(), height=550)
        else:
            st.info("Нажмите 'Построить 3D', чтобы увидеть молекулу")

    with col2:
        st.subheader("Свойства")
        data = get_pubchem_data(smiles)
        if data:
            st.metric("MW", data['mw'])
            st.metric("LogP", data['logp'])
            st.write(f"**Формула:** {data['formula']}")

with tab2:
    st.header("📊 Прогнозирование ADMET свойств")
    
    # Информационный блок
    st.info("""
    **ADMET** (Absorption, Distribution, Metabolism, Excretion, Toxicity) — это комплекс параметров, 
    определяющих судьбу лекарства в организме. На этом этапе мы оцениваем, станет ли молекула реальным лекарством.
    """)

    col_link, col_inst = st.columns([1, 2])

    with col_link:
        st.markdown("### 1. Перейдите в сервис")
        # Кнопка-ссылка на ADMETlab 3.0
        st.link_button("Открыть ADMETlab 3.0 ↗️", "https://admetmesh.scbdd.com/service/evaluation/index")
        st.caption("Рекомендуется использовать браузер Chrome или Edge")

    with col_inst:
        st.markdown("### 2. Краткая инструкция")
        st.markdown("""
        1.  Скопируйте **SMILES** вашей молекулы из первой вкладки.
        2.  Вставьте его в поле ввода на сайте ADMETlab.
        3.  Нажмите кнопку **'Submit'**.
        4.  После завершения расчета найдите кнопку **'Download'** (обычно в формате CSV).
        5.  Вернитесь сюда и загрузите полученный файл ниже.
        """)

    st.divider()

    # Блок загрузки и будущего анализа
    st.markdown("### 3. Загрузка и интерпретация результатов")
    uploaded_file = st.file_uploader("Выберите скачанный CSV файл", type="csv")

    if uploaded_file:
        import pandas as pd
        try:
            df = pd.read_csv(uploaded_file)
            st.success("Файл успешно считан!")

            # Выбираем важные параметры для отображения (названия колонок в ADMETlab 3.0)
            # Примечание: ADMETlab может менять заголовки, поэтому мы ищем совпадения
            
            # 1. Проверка Правила Липинского (Drug-likeness)
            st.subheader("🎯 Оценка Drug-likeness")
            
            # Пытаемся достать данные из CSV (зависит от версии ADMETlab)
            # В ADMETlab 3.0 обычно колонки называются 'MW', 'LogP', 'HBD', 'HBA'
            mw = df['MW'].iloc[0] if 'MW' in df.columns else data.get('mw', 0)
            logp = df['LogP'].iloc[0] if 'LogP' in df.columns else data.get('logp', 0)
            hbd = df['HBD'].iloc[0] if 'HBD' in df.columns else 0
            hba = df['HBA'].iloc[0] if 'HBA' in df.columns else 0

            # Отрисовка индикаторов Липинского
            l1, l2, l3, l4 = st.columns(4)
            l1.metric("MW < 500", f"{mw:.1f}", delta=None if mw <= 500 else "FAIL", delta_color="inverse")
            l2.metric("LogP < 5", f"{logp:.2f}", delta=None if logp <= 5 else "FAIL", delta_color="inverse")
            l3.metric("HBD < 5", int(hbd), delta=None if hbd <= 5 else "FAIL", delta_color="inverse")
            l4.metric("HBA < 10", int(hba), delta=None if hba <= 10 else "FAIL", delta_color="inverse")

            st.divider()

            # 2. Фармакокинетика и Токсичность
            st.subheader("💊 Прогноз Фармакокинетики")
            
            c1, c2 = st.columns(2)
            
            with c1:
                st.write("**Абсорбция и Распределение**")
                # Пример отображения BBB (Blood-Brain Barrier)
                if 'BBB' in df.columns:
                    bbb_val = df['BBB'].iloc[0]
                    color = "green" if bbb_val < 0.5 else "orange"
                    st.markdown(f"Проникновение через ГЭБ (BBB): <span style='color:{color}'>{bbb_val}</span>", unsafe_allow_html=True)
                
                if 'HIA' in df.columns:
                    hia_val = df['HIA'].iloc[0]
                    st.write(f"Всасываемость в кишечнике (HIA): {hia_val}")

            with c2:
                st.write("**Токсичность**")
                # Пример отображения hERG (Кардиотоксичность)
                if 'hERG' in df.columns:
                    herg_val = df['hERG'].iloc[0]
                    risk = "Высокий" if herg_val > 0.7 else "Низкий"
                    st.warning(f"Риск блокировки hERG каналов: {risk} ({herg_val})")

            # 3. Полная таблица для продвинутых
            with st.expander("Открыть полную таблицу результатов"):
                st.dataframe(df.T) # Транспонируем для удобства чтения на смартфонах

        except Exception as e:
            st.error(f"Не удалось проанализировать CSV. Убедитесь, что это файл из ADMETlab. Ошибка: {e}")

with tab3:
    st.header("Образовательный блок")
    st.video("https://www.youtube.com/watch?v=dQw4w9WgXcQ")
    st.markdown("### Контрольный вопрос\nКак влияет LogP на ГЭБ?")
    st.text_area("Ваш ответ")