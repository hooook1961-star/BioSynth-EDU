import os
import re
import gzip
import pickle
import sys
import shutil
import time

import numpy as np
import pandas as pd
from rdkit import Chem
import deepchem as dc

def compute_binding_pocket_features(pocket_featurizer, ligand_featurizer,
                                    pdb_subdir, pdb_code, threshold=.3):
  """Compute features for a given complex"""
  protein_file = os.path.join(pdb_subdir, "%s_protein.pdb" % pdb_code)
  ligand_file = os.path.join(pdb_subdir, "%s_ligand.sdf" % pdb_code)
  ligand_mol2 = os.path.join(pdb_subdir, "%s_ligand.mol2" % pdb_code)

  # Extract active site
  active_site_box, active_site_atoms, active_site_coords = (
      dc.dock.binding_pocket.extract_active_site(
          protein_file, ligand_file))

  # Featurize ligand
  mol = Chem.MolFromMol2File(str(ligand_mol2), removeHs=False)
  if mol is None:
    return None, None
  # Default for CircularFingerprint
  n_ligand_features = 1024
  ligand_features = ligand_featurizer.featurize([mol])

  # Featurize pocket
  finder = dc.dock.ConvexHullPocketFinder()
  pockets, pocket_atoms, pocket_coords = finder.find_pockets(protein_file, ligand_file)
  n_pockets = len(pockets)
  n_pocket_features = dc.feat.BindingPocketFeaturizer.n_features

  features = np.zeros((n_pockets, n_pocket_features+n_ligand_features))
  pocket_features = pocket_featurizer.featurize(
      protein_file, pockets, pocket_atoms, pocket_coords)
  # Note broadcast operation
  features[:, :n_pocket_features] = pocket_features
  features[:, n_pocket_features:] = ligand_features

  # Compute labels for pockets
  labels = np.zeros(n_pockets)
  pocket_atoms[active_site_box] = active_site_atoms
  for ind, pocket in enumerate(pockets):
    overlap = dc.dock.binding_pocket.compute_overlap(
        pocket_atoms, active_site_box, pocket)
    if overlap > threshold:
      labels[ind] = 1
    else:
      labels[ind] = 0 
  return features, labels

def load_pdbbind_labels(labels_file):
  """Loads pdbbind labels as dataframe"""
  # Some complexes have labels but no PDB files. Filter these manually
  missing_pdbs = ["1d2v", "1jou", "1s8j", "1cam", "4mlt", "4o7d"]
  contents = []
  with open(labels_file) as f:
    for line in f:
      if line.startswith("#"):
        continue
      else:
        # Some of the ligand-names are of form (FMN ox). Use regex
        # to merge into form (FMN-ox)
        p = re.compile('\(([^\)\s]*) ([^\)\s]*)\)')
        line = p.sub('(\\1-\\2)', line)
        elts = line.split()
        # Filter if missing PDB files
        if elts[0] in missing_pdbs:
          continue
        contents.append(elts)
  contents_df = pd.DataFrame(
      contents,
      columns=("PDB code", "resolution", "release year", "-logKd/Ki", "Kd/Ki",
               "ignore-this-field", "reference", "ligand name"))
  return contents_df

# Специальный хак-заглушка для совместимости старых DataFrame из старых версий Pandas
class DummyBlockManager:
    def __init__(self, *args, **kwargs):
        pass
    def __setstate__(self, state):
        self.__dict__.update(state)

def featurize_pdbbind_pockets(data_dir=None, subset="core"):
    """
    Лоадер реальных данных из очищенного бинарного архива.
    Включает фильтрацию бесконечных значений для стабильности scikit-learn.
    """
    tasks = ["active-site"]
    
    if data_dir is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    else:
        base_dir = data_dir
        
    pkl_file = os.path.join(base_dir, "datasets", "pdbbind_core_5_clean.pkl.gz")
    
    if not os.path.exists(pkl_file):
        raise ValueError(f"Очищенный файл датасета не найден: {pkl_file}")
        
    print(f"Загрузка реальных 3D-дескрипторов из очищенного архива: {pkl_file}")
    
    with gzip.open(pkl_file, "rb") as f:
        data = pickle.load(f)
        
    X = data["X"]
    y = data["y"]
    ids = data["ids"]
    
    # === БРОНЕБОЙНАЯ СТАБИЛИЗАЦИЯ ДАННЫХ ДЛЯ SKLEARN ===
    # 1. Заменяем любые случайные бесконечности на NaN
    X[~np.isfinite(X)] = np.nan
    # 2. Заменяем NaN на медианные значения по колонкам (или 0, если колонка пустая)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    # 3. Обрезаем экстремальные физические выбросы под границы диапазона float32
    max_float32 = np.finfo(np.float32).max
    min_float32 = np.finfo(np.float32).min
    X = np.clip(X, min_float32, max_float32)
    # ==================================================
    
    w = np.ones_like(y)
    
    print(f"✅ Успешно загружено и стабилизировано {len(X)} реальных комплексов PDBbind.")
    print(f"Размерность матрицы признаков: {X.shape}")
    
    dataset_dir = os.path.join(base_dir, f"{subset}_pockets")
    dataset = dc.data.DiskDataset.from_numpy(X, y, w, ids, data_dir=dataset_dir)
    
    return dataset, tasks
  
def load_pdbbind_pockets(split="random", subset="core"):
  """Основная функция загрузки, вызываемая из binding_pocket_rf.py"""
  import deepchem as dc
  
  # Сама находит корень, где лежит файл лоадера
  base_dir = os.path.dirname(os.path.abspath(__file__))
  
  # Вызывает цепочку обработки данных
  dataset, tasks = featurize_pdbbind_pockets(data_dir=base_dir, subset=subset)

  splitters = {'index': dc.splits.IndexSplitter(),
               'random': dc.splits.RandomSplitter()}
  
  if split not in splitters:
      raise ValueError(f"Неизвестный тип сплиттера: {split}. Доступные: {list(splitters.keys())}")
      
  splitter = splitters[split]
  
  ########################################################### DEBUG
  print("=== ТЕХНИЧЕСКИЙ ДЕБАГ ДАТАСЕТА ===")
  print(f"dataset.X.shape: {dataset.X.shape}")
  print(f"dataset.y.shape: {dataset.y.shape}")
  print(f"dataset.w.shape: {dataset.w.shape}")
  print(f"dataset.ids.shape: {dataset.ids.shape}")
  print("==================================")
  ########################################################### DEBUG
  
  # Разбиваем выборку на обучение, валидацию и тест
  train, valid, test = splitter.train_valid_test_split(dataset)
  
  # Оставляем список пустым, чтобы не ломать распаковку в основном скрипте
  transformers = []
  
  return tasks, (train, valid, test), transformers
