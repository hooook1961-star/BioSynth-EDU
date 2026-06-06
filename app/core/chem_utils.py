import numpy as np
import pubchempy as pcp
import pandas as pd
from functools import lru_cache
from pathlib import Path
import math
import joblib
import pickle
import os
import io
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdFingerprintGenerator
from meeko import MoleculePreparation
from chembl_webresource_client.new_client import new_client

def safe_float(value, default=0.0):
    """
    Безопасно преобразует значение в float. 
    Если значение — строка с запятой, None или ошибка, возвращает default.
    """
    if pd.isna(value): 
        return default
    try:
        clean_val = str(value).replace(',', '.').strip()
        return float(clean_val)
    except (ValueError, TypeError):
        return default

def smiles_to_3d_block(smiles: str, optimize=False) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if optimize:
            AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToMolBlock(mol)
    except: return None

def get_pubchem_data(smiles: str):
    """
    Получает расширенные данные из PubChem.
    """
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if not compounds:
            return None
        
        cmp = compounds[0]
        return {
            "mw": cmp.molecular_weight,
            "logp": cmp.xlogp if cmp.xlogp else "Н/Д",
            "formula": cmp.molecular_formula,
            "iupac": cmp.iupac_name,
            "inchikey": cmp.inchikey,
            "canonical_smiles": cmp.canonical_smiles,
            "rotatable_bonds": cmp.rotatable_bond_count,
            "charge": cmp.charge
        }
    except Exception as e:
        print(f"Ошибка PubChem: {e}")
        return None

def get_chembl_data(inchikey: str):
    """
    Получает данные о механизме действия и статусе лекарства из ChEMBL.
    """
    try:
        molecule = new_client.molecule
        res = molecule.filter(molecule_structures__standard_inchi_key=inchikey).only(['molecule_chembl_id', 'max_phase', 'pref_name', 'indication_class'])
        
        if res:
            mol_data = res[0]
            mechanism = new_client.mechanism
            mechanisms = mechanism.filter(molecule_chembl_id=mol_data['molecule_chembl_id'])
            
            mech_list = []
            for m in mechanisms:
                target = new_client.target.filter(target_chembl_id=m['target_chembl_id']).only(['pref_name'])
                target_name = target[0]['pref_name'] if target else "Unknown Target"
                mech_list.append(f"{m['action_type']} ({target_name})")
            
            return {
                "chembl_id": mol_data['molecule_chembl_id'],
                "pref_name": mol_data.get('pref_name', "N/A"),
                "max_phase": mol_data.get('max_phase', 0),
                "mechanisms": mech_list if mech_list else ["Данные о механизме отсутствуют"]
            }
    except Exception as e:
        print(f"ChEMBL Error: {e}")
    return None
    
def prepare_ligand_for_docking(smiles: str):
    """
    Профессиональная подготовка лиганда в формат PDBQT
    с расчетом торсионов и зарядов.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        mol = Chem.AddHs(mol)
        
        # 1. Генерация 3D и минимизация
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        
        # 2. Подготовка через Meeko
        preparer = MoleculePreparation()
        preparer.prepare(mol)
        
        # 3. Получение PDBQT текста
        pdbqt_block = preparer.write_pdbqt_string()
        return pdbqt_block
    except Exception as e:
        print(f"Meeko Error: {e}")
        return None

@st.cache_resource
def load_scpdb_database():
    file_path = os.path.join("data", "target_database.pkl")
        
    try:
        with open(file_path, "rb") as f:
            return pickle.load(f)
    except Exception as e:
        st.error(f"Критическая ошибка: не удалось загрузить файл базы scPDB по пути {file_path}. Ошибка: {e}")
        return {}

MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def classify_tanimoto_similarity(similarity: float) -> str:
    if similarity >= 0.70:
        return "high"
    if similarity >= 0.50:
        return "moderate"
    if similarity >= 0.35:
        return "weak"
    return "low"


def extract_pdb_id_from_lig_key(custom_name: str) -> str:
    if custom_name.startswith("lig_"):
        return custom_name.removeprefix("lig_").upper()
    return custom_name.upper()


def _clamp(x, lo=0.0, hi=1.0):
    try:
        x = float(x)
    except Exception:
        return 0.0

    if not np.isfinite(x):
        return 0.0

    return max(lo, min(hi, x))


def _safe_float(value, default=0.0):
    try:
        x = float(value)
    except Exception:
        return default

    if not np.isfinite(x):
        return default

    return x


def _exp_similarity(a, b, scale):
    try:
        a = float(a)
        b = float(b)
        scale = float(scale)
    except Exception:
        return 0.0

    if not np.isfinite(a) or not np.isfinite(b) or scale <= 0:
        return 0.0

    return math.exp(-abs(a - b) / scale)


def _get_runtime_index_path() -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    return repo_root / "data" / "scpdb_runtime_index.joblib"


@lru_cache(maxsize=1)
def _load_scpdb_runtime_index():
    path = _get_runtime_index_path()

    if not path.exists():
        return None

    try:
        runtime = joblib.load(path)
    except Exception:
        return None

    required_keys = {"target_keys", "fingerprints", "metadata"}

    if not required_keys.issubset(set(runtime.keys())):
        return None

    return runtime


def _query_descriptors(mol):
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "heavy_atom_count": Descriptors.HeavyAtomCount(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
    }


def _descriptor_compatibility(query_desc, row):
    scores = [
        _exp_similarity(query_desc["mw"], row.get("ligand_mw"), 220.0),
        _exp_similarity(query_desc["logp"], row.get("ligand_logp"), 2.5),
        _exp_similarity(query_desc["tpsa"], row.get("ligand_tpsa"), 90.0),
        _exp_similarity(query_desc["hbd"], row.get("ligand_hbd"), 3.0),
        _exp_similarity(query_desc["hba"], row.get("ligand_hba"), 5.0),
    ]

    ref_charge = _safe_float(row.get("ligand_formal_charge"), 0.0)
    charge_score = 1.0 if query_desc["formal_charge"] == ref_charge else 0.5
    scores.append(charge_score)

    return float(np.mean(scores))


def _site_compatibility(query_desc, row):
    query_polarity = _clamp(query_desc["tpsa"] / 140.0)
    query_hydrophobicity = _clamp(query_desc["logp"] / 5.0)

    site_polar = _safe_float(row.get("site_polar_atom_fraction"), 0.0)
    site_hydro = _safe_float(row.get("site_hydrophobic_atom_fraction"), 0.0)
    cavity6_count = _safe_float(row.get("cavity6_atom_count"), 0.0)

    polar_score = 1.0 - _clamp(abs(query_polarity - site_polar) * 2.0)
    hydro_score = 1.0 - _clamp(abs(query_hydrophobicity - site_hydro) * 2.0)

    required_cavity_points = max(query_desc["heavy_atom_count"] * 3.0, 1.0)
    size_score = _clamp(cavity6_count / required_cavity_points)

    metal_count = _safe_float(row.get("site_metal_atom_count"), 0.0)
    metal_score = 0.6 if metal_count > 0 and query_desc["hba"] < 2 else 1.0

    return float(np.mean([polar_score, hydro_score, size_score, metal_score]))


def _interaction_compatibility(query_desc, row):
    ifp_active_bits = _safe_float(row.get("ifp_active_bits"), 0.0)
    ifp_residues = _safe_float(row.get("ifp_residues_with_interactions"), 0.0)

    ints_sel = _safe_float(row.get("ints_group_SEL"), 0.0)
    ints_all = _safe_float(row.get("ints_group_ALL"), 0.0)
    ints_gll = _safe_float(row.get("ints_group_GLL"), 0.0)
    ints_sep = _safe_float(row.get("ints_group_SEP"), 0.0)
    ints_alp = _safe_float(row.get("ints_group_ALP"), 0.0)
    ints_glp = _safe_float(row.get("ints_group_GLP"), 0.0)
    ints_sec = _safe_float(row.get("ints_group_SEC"), 0.0)
    ints_alc = _safe_float(row.get("ints_group_ALC"), 0.0)
    ints_glc = _safe_float(row.get("ints_group_GLC"), 0.0)

    ifp_score = _clamp(ifp_active_bits / 24.0)
    residue_score = _clamp(ifp_residues / 12.0)

    query_hbond_capacity = query_desc["hbd"] + query_desc["hba"]

    polar_interaction_groups = (
        ints_sel
        + ints_all
        + ints_sep
        + ints_alp
        + ints_sec
        + ints_alc
    )

    hbond_score = _clamp(
        (query_hbond_capacity + 1.0) / (polar_interaction_groups + 1.0)
    )

    hydrophobic_groups = ints_gll + ints_glp + ints_glc
    query_lipophilicity = _clamp(query_desc["logp"] / 5.0)

    hydrophobic_score = _clamp(
        (query_lipophilicity + 0.25) * (1.0 if hydrophobic_groups > 0 else 0.6)
    )

    return float(np.mean([ifp_score, residue_score, hbond_score, hydrophobic_score]))


def _calculate_hybrid_score(
    ligand_similarity,
    descriptor_score,
    site_score,
    interaction_score,
):
    compatibility_score = (
        0.40 * descriptor_score
        + 0.35 * site_score
        + 0.25 * interaction_score
    )

    ligand_gate = _clamp(ligand_similarity / 0.35)

    final_score = (
        0.75 * ligand_similarity
        + 0.25 * compatibility_score * ligand_gate
    )

    return {
        "final_score": float(final_score),
        "compatibility_score": float(compatibility_score),
        "ligand_gate": float(ligand_gate),
    }


def _confidence_from_scores(final_score, ligand_similarity):
    if ligand_similarity < 0.30:
        return "low"

    if final_score >= 0.70 and ligand_similarity >= 0.55:
        return "high"

    if final_score >= 0.55 and ligand_similarity >= 0.40:
        return "moderate"

    if final_score >= 0.40:
        return "weak"

    return "low"


def _format_vina_config_txt(
    receptor_filename: str,
    ligand_filename: str,
    center,
    size,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: int = 3,
) -> str:
    return (
        f"receptor = {receptor_filename}\n"
        f"ligand = {ligand_filename}\n\n"
        f"center_x = {center['x']:.3f}\n"
        f"center_y = {center['y']:.3f}\n"
        f"center_z = {center['z']:.3f}\n\n"
        f"size_x = {size['x']:.3f}\n"
        f"size_y = {size['y']:.3f}\n"
        f"size_z = {size['z']:.3f}\n\n"
        f"exhaustiveness = {int(exhaustiveness)}\n"
        f"num_modes = {int(num_modes)}\n"
        f"energy_range = {int(energy_range)}\n"
    )


def _format_vina_config_yml(
    receptor_filename: str,
    ligand_filename: str,
    center,
    size,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: int = 3,
) -> str:
    return (
        f"receptor: {receptor_filename}\n"
        f"ligand: {ligand_filename}\n"
        f"center:\n"
        f"  x: {center['x']:.3f}\n"
        f"  y: {center['y']:.3f}\n"
        f"  z: {center['z']:.3f}\n"
        f"size:\n"
        f"  x: {size['x']:.3f}\n"
        f"  y: {size['y']:.3f}\n"
        f"  z: {size['z']:.3f}\n"
        f"exhaustiveness: {int(exhaustiveness)}\n"
        f"num_modes: {int(num_modes)}\n"
        f"energy_range: {int(energy_range)}\n"
    )


def _build_docking_payload(row):
    required_columns = [
        "box_center_x",
        "box_center_y",
        "box_center_z",
        "box_size_x",
        "box_size_y",
        "box_size_z",
    ]

    missing_columns = [
        col for col in required_columns
        if col not in row.index
    ]

    if missing_columns:
        return {
            "box_center": None,
            "box_size": None,
            "box_center_source": None,
            "box_size_source": None,
            "docking_box_status": "unavailable",
            "docking_config_txt": "",
            "docking_config_yml": "",
            "receptor_available": False,
            "receptor_source": "external_or_local_required",
            "receptor_filename": "receptor.pdbqt",
            "ligand_filename": "ligand.pdbqt",
        }

    center = {
        "x": _safe_float(row.get("box_center_x"), None),
        "y": _safe_float(row.get("box_center_y"), None),
        "z": _safe_float(row.get("box_center_z"), None),
    }

    size = {
        "x": _safe_float(row.get("box_size_x"), None),
        "y": _safe_float(row.get("box_size_y"), None),
        "z": _safe_float(row.get("box_size_z"), None),
    }

    has_invalid_center = any(v is None for v in center.values())
    has_invalid_size = any(v is None for v in size.values())

    receptor_filename = "receptor.pdbqt"
    ligand_filename = "ligand.pdbqt"

    if has_invalid_center or has_invalid_size:
        return {
            "box_center": None,
            "box_size": None,
            "box_center_source": row.get("box_center_source", None),
            "box_size_source": row.get("box_size_source", None),
            "docking_box_status": "unavailable",
            "docking_config_txt": "",
            "docking_config_yml": "",
            "receptor_available": False,
            "receptor_source": "external_or_local_required",
            "receptor_filename": receptor_filename,
            "ligand_filename": ligand_filename,
        }

    return {
        "box_center": center,
        "box_size": size,
        "box_center_source": row.get("box_center_source", "site.mol2"),
        "box_size_source": row.get("box_size_source", "cavity6.mol2"),
        "docking_box_status": row.get("docking_box_status", "ok"),
        "docking_config_txt": _format_vina_config_txt(
            receptor_filename=receptor_filename,
            ligand_filename=ligand_filename,
            center=center,
            size=size,
        ),
        "docking_config_yml": _format_vina_config_yml(
            receptor_filename=receptor_filename,
            ligand_filename=ligand_filename,
            center=center,
            size=size,
        ),
        "receptor_available": False,
        "receptor_source": "external_or_local_required",
        "receptor_filename": receptor_filename,
        "ligand_filename": ligand_filename,
    }


def _run_legacy_ligand_similarity_screening(
    mol,
    top_n: int = 15,
    min_similarity: float = 0.30,
):
    desc = {}

    target_db = load_scpdb_database()

    if not target_db:
        return {
            "success": False,
            "error_key": "target_error_db_unavailable",
            "desc": desc,
        }

    query_fp = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=2,
        nBits=2048,
    )

    names = list(target_db.keys())
    reference_fps = list(target_db.values())

    try:
        similarities = DataStructs.BulkTanimotoSimilarity(
            query_fp,
            reference_fps,
        )
    except TypeError:
        similarities = [
            DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            for ref_fp in reference_fps
        ]

    candidates = []

    for custom_name, similarity in zip(names, similarities):
        similarity = float(similarity)

        if similarity < min_similarity:
            continue

        pdb_id = extract_pdb_id_from_lig_key(custom_name)
        similarity_level = classify_tanimoto_similarity(similarity)

        candidates.append({
            "id": pdb_id,
            "pdb_id": pdb_id,
            "source_key": custom_name,

            "similarity": similarity,
            "sim": similarity,
            "score": similarity,
            "score_0_100": round(similarity * 100.0, 1),

            "final_score": similarity,
            "confidence": similarity_level,
            "ligand_similarity": similarity,
            "descriptor_compatibility": None,
            "site_compatibility": None,
            "interaction_compatibility": None,

            "box_center": None,
            "box_size": None,
            "box_center_source": None,
            "box_size_source": None,
            "docking_box_status": "legacy_unavailable",
            "docking_config_txt": "",
            "docking_config_yml": "",
            "receptor_available": False,
            "receptor_source": "external_or_local_required",
            "receptor_filename": "receptor.pdbqt",
            "ligand_filename": "ligand.pdbqt",

            "similarity_level": similarity_level,
            "similarity_label_key": f"target_similarity_{similarity_level}",
            "name_key": "target_candidate_name",
            "reason_key": "target_reason_ligand_similarity",
            "interpretation_key": "target_student_interpretation",
            "limitation_key": "target_limitation",
            "method_key": "target_method_short",
        })

    candidates.sort(key=lambda x: x["similarity"], reverse=True)

    if not candidates:
        return {
            "success": True,
            "desc": desc,
            "method": "scPDB ligand similarity screening",
            "method_note_key": "target_method_note",
            "message_key": "target_no_hits_message",
            "raw_score": 0.0,
            "features_expected": 2048,
            "fp_sum": int(query_fp.GetNumOnBits()),
            "top_match": None,
            "all_candidates": [],
            "n_database_entries": len(target_db),
            "n_hits_above_threshold": 0,
            "min_similarity": min_similarity,
            "screening_mode": "legacy_ligand_similarity",
        }

    best_match = candidates[0]

    return {
        "success": True,
        "desc": desc,
        "method": "scPDB ligand similarity screening",
        "method_note_key": "target_method_note",
        "raw_score": best_match["similarity"],
        "features_expected": 2048,
        "fp_sum": int(query_fp.GetNumOnBits()),
        "top_match": best_match,
        "all_candidates": candidates[:top_n],
        "n_database_entries": len(target_db),
        "n_hits_above_threshold": len(candidates),
        "min_similarity": min_similarity,
        "screening_mode": "legacy_ligand_similarity",
    }


def _run_hybrid_scpdb_screening(
    mol,
    top_n: int = 15,
    min_similarity: float = 0.30,
    candidate_pool: int = 1000,
):
    runtime = _load_scpdb_runtime_index()

    if runtime is None:
        return _run_legacy_ligand_similarity_screening(
            mol=mol,
            top_n=top_n,
            min_similarity=min_similarity,
        )

    query_desc = _query_descriptors(mol)
    query_fp = MORGAN_GEN.GetFingerprint(mol)

    metadata = runtime["metadata"].copy()

    similarities = DataStructs.BulkTanimotoSimilarity(
        query_fp,
        runtime["fingerprints"],
    )

    metadata["ligand_similarity"] = similarities

    n_hits_above_threshold = int(
        (metadata["ligand_similarity"] >= min_similarity).sum()
    )

    candidates = (
        metadata
        .sort_values("ligand_similarity", ascending=False)
        .head(candidate_pool)
        .copy()
    )

    rows = []

    for _, row in candidates.iterrows():
        ligand_similarity = _safe_float(row.get("ligand_similarity"), 0.0)

        descriptor_score = _descriptor_compatibility(query_desc, row)
        site_score = _site_compatibility(query_desc, row)
        interaction_score = _interaction_compatibility(query_desc, row)

        score_info = _calculate_hybrid_score(
            ligand_similarity=ligand_similarity,
            descriptor_score=descriptor_score,
            site_score=site_score,
            interaction_score=interaction_score,
        )

        final_score = score_info["final_score"]
        confidence = _confidence_from_scores(final_score, ligand_similarity)
        similarity_level = classify_tanimoto_similarity(ligand_similarity)

        pdb_id = str(row["pdb_id"]).upper()
        entry_id = str(row["entry_id"])
        target_key = str(row["target_key"])

        docking_payload = _build_docking_payload(row)

        rows.append({
            "id": pdb_id,
            "pdb_id": pdb_id,
            "entry_id": entry_id,
            "target_key": target_key,
            "source_key": target_key,

            "similarity": ligand_similarity,
            "sim": ligand_similarity,
            "ligand_similarity": ligand_similarity,

            "score": final_score,
            "score_0_100": round(final_score * 100.0, 1),
            "final_score": final_score,
            "confidence": confidence,

            "descriptor_compatibility": descriptor_score,
            "site_compatibility": site_score,
            "interaction_compatibility": interaction_score,
            "compatibility_score": score_info["compatibility_score"],
            "ligand_gate": score_info["ligand_gate"],

            "similarity_level": similarity_level,
            "similarity_label_key": f"target_similarity_{similarity_level}",
            "name_key": "target_candidate_name",
            "reason_key": "target_reason_ligand_similarity",
            "interpretation_key": "target_student_interpretation",
            "limitation_key": "target_limitation",
            "method_key": "target_method_short",

            **docking_payload,
        })

    rows.sort(key=lambda x: x["final_score"], reverse=True)

    top_candidates = rows[:top_n]
    top_match = top_candidates[0] if top_candidates else None

    if top_match is None:
        return {
            "success": True,
            "desc": {},
            "method": "hybrid_scpdb_target_hypothesis_screening",
            "method_note_key": "target_method_note",
            "message_key": "target_no_hits_message",
            "raw_score": 0.0,
            "features_expected": 2048,
            "fp_sum": int(query_fp.GetNumOnBits()),
            "top_match": None,
            "all_candidates": [],
            "n_database_entries": int(len(metadata)),
            "n_hits_above_threshold": n_hits_above_threshold,
            "min_similarity": min_similarity,
            "screening_mode": "hybrid_scpdb",
        }

    return {
        "success": True,
        "desc": query_desc,
        "method": "hybrid_scpdb_target_hypothesis_screening",
        "method_note_key": "target_method_note",

        # main.py показывает raw_score как лучшее Tanimoto-сходство.
        "raw_score": top_match["similarity"],

        "features_expected": 2048,
        "fp_sum": int(query_fp.GetNumOnBits()),
        "top_match": top_match,
        "all_candidates": top_candidates,
        "n_database_entries": int(len(metadata)),
        "n_hits_above_threshold": n_hits_above_threshold,
        "min_similarity": min_similarity,
        "screening_mode": "hybrid_scpdb",
    }


def run_ai_target_screening(
    smiles_str,
    top_n: int = 15,
    min_similarity: float = 0.30,
):
    desc = {}

    mol = Chem.MolFromSmiles(smiles_str.strip()) if smiles_str else None

    if mol is None:
        return {
            "success": False,
            "error_key": "target_error_invalid_smiles",
            "desc": desc,
        }

    return _run_hybrid_scpdb_screening(
        mol=mol,
        top_n=top_n,
        min_similarity=min_similarity,
        candidate_pool=1000,
    )

def get_quantum_descriptors(smiles: str) -> pd.DataFrame:
    """
    Расчет фундаментальных структурно-химических и полярных дескрипторов.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        Chem.SanitizeMol(mol)
        
        tpsa = Descriptors.TPSA(mol)
        labute_asa = Descriptors.LabuteASA(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        rot_bonds = Descriptors.RotatableBonds(mol)
        steric_index = round(rot_bonds / heavy_atoms, 3) if heavy_atoms > 0 else 0
        hetero_atoms = len([at for at in mol.GetAtoms() if at.GetSymbol() in ['N', 'O']])

        df_data = {
            "Молекулярный дескриптор": [
                "Топологическая полярная поверхность (TPSA, Å²)",
                "Доступный объем Ван-дер-Ваальса (Labute ASA)",
                "Стерический индекс жесткости (RotB/Heavy)",
                "Количество гетероатомов-доноров (N, O)"
            ],
            "Значение": [
                f"{tpsa:.2f} Å²",
                f"{labute_asa:.2f}",
                f"{steric_index}",
                f"{hetero_atoms}"
            ]
        }
        return pd.DataFrame(df_data)
    except Exception as e:
        print(f"Ошибка дескрипторов: {e}")
        return None
