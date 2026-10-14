"""Offline unit tests: execute the unchanged Vina helpers without UI/network imports.

Run: python -m unittest discover -s tests -p test_pdb_search.py
Requires numpy and pandas, already listed by the application.
"""
import ast
import math
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

APP = Path(__file__).resolve().parents[1] / 'app'
sys.path.insert(0, str(APP))
from core.pdb_search import PDBSearchError, normalize_pdb_id, search_protein_by_pdb_id


def existing_helpers():
    source = ast.parse((APP / 'core/chem_utils.py').read_text(encoding='utf-8'))
    names = {'_build_docking_payload', '_safe_float', '_format_vina_config_txt',
             '_format_vina_config_yml', 'extract_pdb_id_from_lig_key'}
    functions = [n for n in source.body if isinstance(n, ast.FunctionDef) and n.name in names]
    module = types.ModuleType('core.chem_utils')
    module.__dict__.update(np=np, math=math)
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(APP / 'core/chem_utils.py'), 'exec'), module.__dict__)
    return module


class PDBSearchTests(unittest.TestCase):
    def setUp(self):
        self.chem = existing_helpers()
        self.rows = [dict(pdb_id='4a5s', entry_id='4a5s_1', target_key='lig_4a5s',
                          box_center_x=1.25, box_center_y=-2, box_center_z=0,
                          box_size_x=20, box_size_y=21, box_size_z=22,
                          box_center_source='site.mol2', box_size_source='cavity6.mol2')]
        self.chem._load_scpdb_runtime_index = lambda: {'metadata': pd.DataFrame(self.rows)}
        self.chem.load_scpdb_database = lambda: {}
        self.modules = patch.dict(sys.modules, {'core.chem_utils': self.chem})
        self.modules.start()
        self.addCleanup(self.modules.stop)

    def test_normalization(self):
        self.assertEqual(normalize_pdb_id(' 4a5s '), '4A5S')
        for value in ['', 'AAAA', '4A5S/../', '4A5', '4A5SS']:
            with self.assertRaises(PDBSearchError):
                normalize_pdb_id(value)

    def test_same_coordinates_and_config_as_existing_mode(self):
        original = self.chem._build_docking_payload(pd.Series(self.rows[0]))
        result = search_protein_by_pdb_id(' 4a5s ')['all_candidates'][0]
        for key, value in original.items():
            self.assertEqual(result[key], value, key)
        self.assertIn('center_x = 1.250', result['docking_config_txt'])
        self.assertIn('size_y = 21.000', result['docking_config_txt'])

    def test_multiple_sites_and_exact_match(self):
        self.rows += [dict(self.rows[0], entry_id='4a5s_2'), dict(self.rows[0], pdb_id='4A5S_extra')]
        self.assertEqual(len(search_protein_by_pdb_id('4A5S')['all_candidates']), 2)

    def test_not_found(self):
        with self.assertRaises(PDBSearchError):
            search_protein_by_pdb_id('1ABC')

    def test_missing_coordinates(self):
        del self.rows[0]['box_center_x']
        self.assertEqual(search_protein_by_pdb_id('4A5S')['all_candidates'][0]['docking_config_txt'], '')

    def test_invalid_coordinates_and_sizes(self):
        for field, value in [('box_center_x', math.nan), ('box_center_y', math.inf),
                             ('box_size_z', 0), ('box_size_x', -5), ('box_size_y', 'invalid')]:
            with self.subTest(field=field, value=value):
                original = self.rows[0][field]
                self.rows[0][field] = value
                self.assertEqual(search_protein_by_pdb_id('4A5S')['all_candidates'][0]['docking_config_txt'], '')
                self.rows[0][field] = original

    def test_legacy_fallback_has_no_fabricated_box(self):
        self.chem._load_scpdb_runtime_index = lambda: None
        self.chem.load_scpdb_database = lambda: {'lig_4a5s': object()}
        result = search_protein_by_pdb_id('4A5S')['all_candidates'][0]
        self.assertIsNone(result['box_center'])
        self.assertEqual(result['docking_config_txt'], '')

    def test_unavailable_database(self):
        self.chem._load_scpdb_runtime_index = lambda: None
        with self.assertRaises(PDBSearchError):
            search_protein_by_pdb_id('4A5S')

    def test_invalid_schema(self):
        del self.rows[0]['pdb_id']
        with self.assertRaises(PDBSearchError):
            search_protein_by_pdb_id('4A5S')


if __name__ == '__main__':
    unittest.main()
