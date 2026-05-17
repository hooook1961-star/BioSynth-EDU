import numpy as np
import pubchempy as pcp
import pandas as pd
import pickle
import os
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
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

def run_ai_target_screening(smiles_str):
    desc = {} 
    
    mol = Chem.MolFromSmiles(smiles_str) if smiles_str else None
    if mol is None:
        return {"error": "Некорректный SMILES", "desc": desc}

    target_db = load_scpdb_database()
    if not target_db:
        return {"error": "База данных мишеней недоступна или пуста", "desc": desc}

    query_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    
    scored_proteins = []

    for custom_name, ref_fp in target_db.items():
        similarity = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
        
        pdb_id = custom_name.replace("lig_", "").upper()
        dynamic_score = 5.0 + (similarity * 4.0)

        scored_proteins.append({
            "id": pdb_id,
            "name": f"Мишень из базы scPDB (ID: {pdb_id})",
            "reason": f"Индекс сходства Танимото с нативным лигандом: {similarity:.2f}.",
            "score": float(dynamic_score),
            "sim": similarity
        })

    scored_proteins.sort(key=lambda x: x["sim"], reverse=True)
    best_match = scored_proteins[0]

    return {
        "success": True,
        "desc": desc,                 
        "raw_score": best_match["sim"], 
        "features_expected": 2048,
        "fp_sum": int(np.array(query_fp).sum()),
        "top_match": best_match,
        "all_candidates": scored_proteins[:15] 
    }
