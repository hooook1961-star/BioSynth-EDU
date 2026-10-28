import hashlib
import logging

import streamlit as st

from core.complex_analysis import (ComplexAnalysisError, KINDS, analyze_complex,
                                   contacts_csv, inspect_pdb, table_rows)
from core.interaction_diagram import render_diagram


def render_complex_analysis():
    st.subheader('Анализ комплекса белок–лиганд')
    st.write('Загрузите готовый PDB-комплекс с одной позой лиганда после AutoDock Vina или SwissDock.')
    uploaded = st.file_uploader('Комплекс в формате PDB', type=['pdb'], key='complex_pdb_upload')
    if uploaded is None:
        st.session_state.pop('complex_analysis_result', None)
        st.session_state.pop('complex_analysis_input', None)
        st.session_state.pop('complex_diagram', None)
        return
    data = uploaded.getvalue()
    digest = hashlib.sha256(data).hexdigest()
    if st.session_state.get('complex_analysis_upload') != digest:
        st.session_state.pop('complex_analysis_result', None)
        st.session_state.pop('complex_analysis_input', None)
        st.session_state['complex_analysis_upload'] = digest
    try:
        inspected = inspect_pdb(data)
    except ComplexAnalysisError as exc:
        st.error(str(exc))
        return
    candidates = inspected['candidates']
    if len(candidates) == 1:
        selected = candidates[0]['key']
        st.caption(f'Лиганд: {selected} · тяжёлых атомов: {candidates[0]["atom_count"]}')
    else:
        st.info('В комплексе несколько органических компонентов. Укажите анализируемый лиганд.')
        selected = st.selectbox('Лиганд комплекса', [None] + [c['key'] for c in candidates],
                                format_func=lambda x: x or 'Выберите компонент', key=f'complex_ligand_{digest}')
        if selected is None:
            return
    input_key = (digest, selected)
    if st.session_state.get('complex_analysis_input') != input_key:
        st.session_state.pop('complex_analysis_result', None)
    if st.button('Анализировать взаимодействия', key='analyze_complex_button'):
        st.session_state.pop('complex_analysis_result', None)
        st.session_state.pop('complex_diagram', None)
        st.session_state['complex_analysis_input'] = input_key
        try:
            with st.spinner('Определение взаимодействий белка и лиганда…'):
                st.session_state['complex_analysis_result'] = analyze_complex(data, selected)
        except ComplexAnalysisError as exc:
            st.error(str(exc))
        except Exception:
            logging.getLogger(__name__).exception('Complex interaction analysis failed')
            st.error('Не удалось завершить анализ комплекса. Проверьте подготовку PDB-файла.')
    result = st.session_state.get('complex_analysis_result')
    if result is None:
        return
    st.caption(f"PLIP {result['plip_version']}. Химические группы и водороды определяются автоматически; проверьте подготовку и протонирование комплекса.")
    rows = table_rows(result)
    if rows:
        st.dataframe(rows, hide_index=True, use_container_width=True)
    else:
        st.info('Водородные, гидрофобные, солевые и ароматические взаимодействия по критериям PLIP не обнаружены.')
    st.download_button('Скачать таблицу CSV', contacts_csv(result), 'interactions.csv',
                       mime='text/csv', key='complex_csv_download')
    chosen = st.multiselect('Типы взаимодействий на диаграмме', list(KINDS),
                            default=list(KINDS), format_func=lambda k: KINDS[k][0], key='complex_diagram_kinds')
    st.caption('Фильтр меняет только рисунок. CSV содержит полную таблицу взаимодействий.')
    try:
        # Session-only cache: no cross-user storage of uploaded structures.
        diagram_key = (input_key, tuple(chosen))
        cached = st.session_state.get('complex_diagram')
        if not cached or cached[0] != diagram_key:
            cached = (diagram_key, render_diagram(result, chosen))
            st.session_state['complex_diagram'] = cached
        svg, png = cached[1]
        st.image(png, caption='2D-диаграмма взаимодействий выбранного лиганда', use_container_width=True)
        left, right = st.columns(2)
        left.download_button('Скачать PNG', png, 'interactions.png', mime='image/png', key='complex_png_download')
        right.download_button('Скачать SVG', svg, 'interactions.svg', mime='image/svg+xml', key='complex_svg_download')
    except ComplexAnalysisError as exc:
        st.warning(str(exc))
    except Exception:
        logging.getLogger(__name__).exception('Interaction diagram failed')
        st.warning('Таблица готова, но построить химически корректную диаграмму не удалось. Проверьте структуру лиганда и зависимости визуализации.')
    st.caption('Взаимодействия характеризуют геометрию выбранной позы и сами по себе не доказывают силу связывания.')
