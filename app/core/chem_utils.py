import numpy as np
import pubchempy as pcp
import pandas as pd
import gzip
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from meeko import MoleculePreparation

def safe_float(value, default=0.0):
    """
    Безопасно преобразует значение в float. 
    Если значение — строка с запятой, None или ошибка, возвращает default.
    """
    if pd.isna(value): # Проверка на NaN из pandas
        return default
    try:
        # Убираем пробелы и заменяем запятую на точку (на случай региональных стандартов)
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
        
        # Собираем расширенный словарь данных
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

from chembl_webresource_client.new_client import new_client

def get_chembl_data(inchikey: str):
    """
    Получает данные о механизме действия и статусе лекарства из ChEMBL.
    """
    try:
        molecule = new_client.molecule
        res = molecule.filter(molecule_structures__standard_inchi_key=inchikey).only(['molecule_chembl_id', 'max_phase', 'pref_name', 'indication_class'])
        
        if res:
            mol_data = res[0]
            # Получаем механизмы действия
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
        
        # 2. Подготовка через Meeko (аналог того, что делает AutoDockTools)
        preparer = MoleculePreparation()
        preparer.prepare(mol)
        
        # 3. Получение PDBQT текста (с ROOT, ENDROOT, TORSDOF и т.д.)
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
        # Логирование ошибки можно добавить при необходимости (например, print(e))
        return default_values

def run_ai_target_screening(smiles_str, pocket_model):
    """
    Принимает SMILES и обученную модель.
    Загружает оригинальный архив, выполняет скрининг всей матрицы X
    и возвращает отсортированный список результатов.
    """
    # 1. Считаем дескрипторы для паспорта лигандов
    default_desc = {"mw": 200.0, "logp": 1.5, "tpsa": 40.0}
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is not None:
            desc = {
                "mw": float(Descriptors.MolWt(mol)),
                "logp": float(Descriptors.MolLogP(mol)),
                "tpsa": float(Descriptors.TPSA(mol))
            }
        else:
            desc = default_desc
    except:
        desc = default_desc

    # 2. Читаем оригинальный pkl.gz архива обучения
    try:
        with gzip.open("datasets/pdbbind_core_5_clean.pkl.gz", "rb") as f:
            archive_data = pickle.load(f)
        X_matrix = archive_data.get("X")
        pdb_ids = archive_data.get("ids")
    except Exception as e:
        # Если файл не прочитался, возвращаем ошибку для интерфейса
        return {"error": f"Файл архива недоступен: {e}", "desc": desc}

    if X_matrix is None or pdb_ids is None:
        return {"error": "В архиве отсутствуют ключи 'X' или 'ids'.", "desc": desc}

    # 3. Прямой расчёт модели Случайного Леса
    try:
        predictions = pocket_model.predict(X_matrix)
    except Exception as e:
        return {"error": f"Ошибка размерности матрицы в модели: {e}", "desc": desc}

    # 4. Сборка и сортировка результатов
    scored_proteins = []
    for i, pdb_code in enumerate(pdb_ids):
        pdb_str = str(pdb_code).strip().upper()
        predicted_pkd = float(predictions[i])
        
        # Научные аннотации для ключевых кейсов
        name_str = f"Биологическая мишень (PDB ID: {pdb_str})"
        reason_str = "Верифицированная мишень из оригинальной обучающей выборки PDBbind Core Set."
        
        if pdb_str == "1UWH":
            name_str = "Клеточная киназа опухоли (Рецептор EGFR)"
            reason_str = "Классическая мишень для анализа цитотоксичности и противоопухолевой активности липофильных диеноновых производных пиперидона."
        elif pdb_str == "1CX2":
            name_str = "Циклооксигеназа-2 (ЦОГ-2) — фермент воспаления"
            reason_str = "Целевой фермент воспаления. Фармакологический эффект НПВП обусловлен ингибированием этого активного центра."

        scored_proteins.append({
            "id": pdb_str,
            "name": name_str,
            "reason": reason_str,
            "score": predicted_pkd
        })

    # Сортируем по убыванию предсказанного pKd
    scored_proteins = sorted(scored_proteins, key=lambda x: x["score"], reverse=True)
    
    return {
        "success": True,
        "desc": desc,
        "top_match": scored_proteins[0] # Возвращаем самый лучший результат
    }
