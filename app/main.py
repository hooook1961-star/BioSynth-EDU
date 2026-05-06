import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import pandas as pd
from core.chem_utils import smiles_to_3d_block, get_pubchem_data, get_chembl_data, prepare_ligand_for_docking

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

# Вкладка ADMET 
with tab2:
    st.header("📊 ADMET Анализ и интерпретация")
    
    # Блок с инструкциями и кнопками
    st.markdown("""
    Для проведения глубокого анализа ADMET:
    1. Нажмите кнопку **«Открыть ADMETlab 3.0»** ниже.
    2. Скопируйте ваш SMILES (из первой вкладки) и вставьте его в поле ввода на сайте.
    3. После завершения расчета скачайте результат в формате **CSV**.
    4. Загрузите файл сюда для автоматической интерпретации.
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
                    if -1 < val < 5: # Расширили диапазон для полярных соединений
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
            1. Добавление неявных водородов (H-atoms).
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
                
                # --- НОВЫЙ БЛОК: Визуализация подготовленного лиганда ---
                st.write("🔍 **Просмотр подготовленной структуры:**")
                import py3Dmol
                
                view = py3Dmol.view(width=300, height=300)
                # Загружаем именно PDBQT данные
                view.addModel(st.session_state.prepared_pdbqt, 'pdbqt')
                view.setStyle({'stick': {'color': 'spectrum', 'radius': 0.2}, 'sphere': {'scale': 0.3}})
                view.zoomTo()
                
                # Рендерим компонент в Streamlit
                st.components.v1.html(view._make_html(), height=310)
                st.caption("Подготовленный лиганд с водородами (H) и оптимизированными связями.")
            
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
    st.header("Курс лекций и тесты")
    st.video("https://www.youtube.com/watch?v=dQw4w9WgXcQ")
    st.markdown("### Мини-тест\nКак изменится липофильность (LogP) при добавлении -OH группы?")
    st.radio("Выберите ответ:", ["Увеличится", "Уменьшится", "Не изменится"])
