import numpy as np
import pubchempy as pcp
import pandas as pd
import gzip
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
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

def calculate_molecule_descriptors(smiles_str):
    """
    Принимает SMILES и рассчитывает ключевые хемоинформатические дескрипторы.
    Возвращает словарь со значениями по умолчанию в случае ошибки или пустого ввода.
    """
    default_values = {"mw": 200.0, "logp": 1.5, "tpsa": 40.0, "hbd": 1, "hba": 2}
    if not smiles_str:
        return default_values
        
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return default_values
            
        return {
            "mw": float(Descriptors.MolWt(mol)),
            "logp": float(Descriptors.MolLogP(mol)),
            "tpsa": float(Descriptors.TPSA(mol)),
            "hbd": int(Descriptors.NumHDonors(mol)),
            "hba": int(Descriptors.NumHAcceptors(mol))
        }
    except Exception as e:
        return default_values

def run_ai_target_screening(smiles_str, pocket_model):
    """
    Абсолютно голый диагностический QSAR без цензуры и заглушек.
    """
    desc = calculate_molecule_descriptors(smiles_str)
    mol = Chem.MolFromSmiles(smiles_str) if smiles_str else None
    if mol is None:
        return {"error": "Некорректный SMILES", "desc": desc}

    # 1. Генерируем фингерпринт
    try:
        required_features = pocket_model.n_features_in_
    except:
        required_features = 2048 

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=required_features)
    fp_array = np.zeros((1,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, fp_array)
    ligand_vector = fp_array.reshape(1, -1).astype(np.float32)

    # 2. Берем сырой прогноз модели
    try:
        raw_prediction = pocket_model.predict(ligand_vector)
        raw_score = float(raw_prediction[0])
    except Exception as e:
        return {"error": f"Ошибка .predict(): {e}", "desc": desc}

    # Возвращаем СЫРЫЕ данные для диагностики
    return {
        "success": True,
        "desc": desc,
        "raw_score": raw_score,
        "features_expected": required_features,
        "fp_sum": int(fp_array.sum()) # Показывает, сколько единичек в фингерпринте
    }
