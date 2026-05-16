"""
PDBBind binding pocket dataset loader.
"""

from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

import os
import numpy as np
import pandas as pd
import shutil
import time
import re
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

def featurize_pdbbind_pockets(data_dir=None, subset="core"):
  """Адаптированный лоадер: жесткие абсолютные пути для Streamlit Cloud"""
  tasks = ["active-site"]
  
  # Если data_dir не передан, берем папку файла, иначе используем переданный абсолютный путь
  if data_dir is None:
      base_dir = os.path.dirname(os.path.abspath(__file__))
  else:
      base_dir = data_dir
  
  # Строим железный путь к архиву
  pkl_file = os.path.join(base_dir, "datasets", "pdbbind_core_5_df.pkl.gz")
  
  if not os.path.exists(pkl_file):
    raise ValueError(f"Файл датасета не найден по пути: {pkl_file}. Проверь папку datasets в корне.")

  print(f"Загрузка и распаковка готового датасета из: {pkl_file}")
  
  try:
      import gzip
      import pickle
      with gzip.open(pkl_file, "rb") as f:
          df = pickle.load(f, encoding="latin1")
  except Exception as e:
      raise RuntimeError(f"Не удалось прочитать .pkl.gz архив. Ошибка: {str(e)}")

  print(f"Успешно загружено {len(df)} complexes.")

  X = np.vstack(df["X"].values)
  y = np.concatenate(df["y"].values)
  w = np.ones_like(y)
  
  all_ids = []
  for _, row in df.iterrows():
      pdb_code = row["pdb_id"]
      p_ids = np.array(["%s%d" % (pdb_code, i) for i in range(len(row["y"]))])
      all_ids.append(p_ids)
  ids = np.concatenate(all_ids)
   
  # Папку кэша DiskDataset создаем тоже строго в абсолютном пути бэкенда
  dataset_dir = os.path.join(base_dir, "%s_pockets" % (subset))
  dataset = dc.data.DiskDataset.from_numpy(X, y, w, ids, data_dir=dataset_dir)
  
  return dataset, tasks

def load_pdbbind_pockets(split="index", subset="core"):
  """Load PDBBind datasets. Does not do train/test split"""
  dataset, tasks = featurize_pdbbind_pockets(subset=subset)

  splitters = {'index': dc.splits.IndexSplitter(),
               'random': dc.splits.RandomSplitter()}
  splitter = splitters[split]
  ########################################################### DEBUG
  print("dataset.X.shape")
  print(dataset.X.shape)
  print("dataset.y.shape")
  print(dataset.y.shape)
  print("dataset.w.shape")
  print(dataset.w.shape)
  print("dataset.ids.shape")
  print(dataset.ids.shape)
  ########################################################### DEBUG
  train, valid, test = splitter.train_valid_test_split(dataset)

  transformers = []
  for transformer in transformers:
    train = transformer.transform(train)
  for transformer in transformers:
    valid = transformer.transform(valid)
  for transformer in transformers:
    test = transformer.transform(test)
  
  return tasks, (train, valid, test), transformers
