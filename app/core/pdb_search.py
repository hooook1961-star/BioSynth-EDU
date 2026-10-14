"""Exact PDB ID lookup using the existing scPDB index and Vina generator."""
import math
import re


class PDBSearchError(ValueError):
    """An actionable lookup error safe to display in the UI."""


def normalize_pdb_id(value):
    pdb_id = str(value).strip().upper()
    if not re.fullmatch(r"[1-9][A-Z0-9]{3}", pdb_id):
        raise PDBSearchError("Введите PDB ID из четырёх символов, например 4A5S.")
    return pdb_id


def search_protein_by_pdb_id(value):
    from core.chem_utils import (
        _load_scpdb_runtime_index, _build_docking_payload,
        load_scpdb_database, extract_pdb_id_from_lig_key,
    )
    pdb_id = normalize_pdb_id(value)
    runtime = _load_scpdb_runtime_index()
    candidates = []
    if runtime is not None:
        metadata = runtime["metadata"]
        if "pdb_id" not in metadata.columns:
            raise PDBSearchError("Индекс scPDB повреждён: отсутствует поле pdb_id.")
        rows = metadata.loc[metadata["pdb_id"].astype(str).str.strip().str.upper() == pdb_id]
        for _, row in rows.iterrows():
            payload = _build_docking_payload(row)
            # Reject invalid database boxes without changing similarity mode.
            center, size = payload.get("box_center"), payload.get("box_size")
            valid = (center and size
                     and all(math.isfinite(center[a]) for a in "xyz")
                     and all(math.isfinite(size[a]) and size[a] > 0 for a in "xyz"))
            if not valid:
                payload.update(box_center=None, box_size=None,
                               docking_config_txt="", docking_config_yml="",
                               docking_box_status="unavailable")
            candidates.append({"pdb_id": pdb_id, "entry_id": str(row.get("entry_id", pdb_id)),
                               "target_key": str(row.get("target_key", "")), **payload})
    else:
        # Same fallback as similarity search. Fingerprints have no coordinates.
        legacy = load_scpdb_database()
        if not legacy:
            raise PDBSearchError("База scPDB недоступна. Проверьте файлы в папке data.")
        for key in legacy:
            if extract_pdb_id_from_lig_key(key).upper() == pdb_id:
                candidates.append({"pdb_id": pdb_id, "entry_id": str(key),
                                   "box_center": None, "box_size": None,
                                   "docking_config_txt": "", "docking_config_yml": "",
                                   "docking_box_status": "legacy_unavailable"})
    if not candidates:
        raise PDBSearchError(f"{pdb_id} не найден в локальной базе scPDB.")
    return {"pdb_id": pdb_id, "all_candidates": candidates}
