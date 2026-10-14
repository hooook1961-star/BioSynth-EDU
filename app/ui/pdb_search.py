import logging

import streamlit as st

from core.pdb_search import PDBSearchError, search_protein_by_pdb_id


def render_pdb_search():
    st.subheader("Поиск белка по PDB ID")
    st.caption("Поиск в той же базе scPDB, что и поиск по сходству. Координаты и размеры бокса берутся из существующей записи базы.")
    value = st.text_input("PDB ID", placeholder="4A5S", key="pdb_id_input")
    if st.button("Найти белок", key="pdb_search_button"):
        st.session_state.pop("pdb_search_result", None)
        try:
            with st.spinner("Поиск в scPDB…"):
                st.session_state.pdb_search_result = search_protein_by_pdb_id(value)
        except PDBSearchError as exc:
            st.error(str(exc))
        except Exception:
            logging.getLogger(__name__).exception("PDB ID search failed")
            st.error("Не удалось прочитать базу scPDB. Проверьте файлы базы и повторите поиск.")
    result = st.session_state.get("pdb_search_result")
    if not result or value.strip().upper() != result["pdb_id"]:
        return
    pdb_id = result["pdb_id"]
    candidates = result["all_candidates"]
    st.success(f"{pdb_id}: найдено записей scPDB — {len(candidates)}.")
    choices = [None] + list(range(len(candidates))) if len(candidates) > 1 else [0]
    selected = st.selectbox("Запись scPDB / сайт связывания", choices,
                            format_func=lambda i: "Выберите запись" if i is None else candidates[i]["entry_id"],
                            key=f"pdb_entry_{pdb_id}")
    if selected is None:
        st.info("Для этого PDB ID найдено несколько записей. Выберите сайт связывания.")
        return
    candidate = candidates[selected]
    config = candidate.get("docking_config_txt", "")
    if not config:
        st.warning("Для этой записи нет корректных координат и размеров бокса. Конфиг не сформирован. В старой базе отпечатков координаты не хранятся.")
        return
    st.caption(f"Источник центра: {candidate.get('box_center_source', '—')}; источник размеров: {candidate.get('box_size_source', '—')}.")
    st.code(config, language="ini")
    st.download_button("Скачать конфиг AutoDock Vina", config, f"{pdb_id}_vina_config.txt", mime="text/plain", key="pdb_config_download")
    st.download_button("Скачать YAML-конфиг", candidate["docking_config_yml"], f"{pdb_id}_vina_config.yml", mime="text/yaml", key="pdb_yaml_download")
    st.info("Конфиг использует receptor.pdbqt и ligand.pdbqt, как в режиме поиска по сходству. Эти файлы нужно подготовить отдельно с сохранением исходной системы координат.")
