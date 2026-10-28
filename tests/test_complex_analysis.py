"""Offline contract tests. Synthetic coordinates are test data, not results."""
import ast
import importlib.util
from pathlib import Path
import subprocess
import sys
import types
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'app'))
from core import complex_analysis as ca
from core.interaction_diagram import diagram_layout


def atom(serial, name, residue, number, element, record='HETATM', chain='A', alt=' ', insertion=' '):
    return (f'{record:<6}{serial:5d} {name:^4}{alt}{residue:>3} {chain}{number:4d}{insertion}   '
            f'{serial * 1.2:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2}  ')


def pdb():
    return ('\n'.join([atom(1, 'CA', 'ALA', 1, 'C', 'ATOM'), atom(2, 'C1', 'LIG', 2, 'C'),
                      atom(3, 'C2', 'LIG', 2, 'C'), atom(4, 'O1', 'LIG', 2, 'O')]) + '\nEND\n').encode()


def xml(kind='hydrogen_bond', extra='<protisdon>True</protisdon><acceptoridx>4</acceptoridx><donoridx>1</donoridx>', distance='dist_d-a', with_contact=True):
    contact = f'<{kind}><restype>ALA</restype><reschain>A</reschain><resnr>1</resnr><{distance}>3.1</{distance}>{extra}</{kind}>' if with_contact else ''
    return f'''<report><plipversion>2.4.0</plipversion><bindingsite>
      <identifiers><hetid>LIG</hetid><chain>A</chain><position>2</position><composite>False</composite><smiles>CCO</smiles></identifiers>
      <interactions><contacts>{contact}</contacts></interactions>
      <mappings><smiles_to_pdb>1:2,2:3,3:4</smiles_to_pdb></mappings>
      </bindingsite></report>'''


class ValidationTests(unittest.TestCase):
    def test_valid_complex_and_water_exclusion(self):
        data = pdb() + (atom(5, 'O', 'HOH', 3, 'O') + '\n').encode()
        self.assertEqual(ca.inspect_pdb(data)['candidates'], [dict(key='LIG:A:2', atom_count=3)])

    def test_multiple_ligands_remain_separate(self):
        data = pdb() + ('\n'.join([atom(5, 'C1', 'DRG', 4, 'C'), atom(6, 'O1', 'DRG', 4, 'O')])).encode()
        self.assertEqual(len(ca.inspect_pdb(data)['candidates']), 2)

    def test_no_protein(self):
        with self.assertRaises(ca.ComplexAnalysisError):
            ca.inspect_pdb((atom(1, 'C1', 'LIG', 2, 'C') + '\n').encode())

    def test_missing_ligand_and_additives(self):
        for name in ['HOH', 'GOL', 'SO4', 'ZN', 'EDO']:
            data = pdb().replace(b'LIG', name.encode().rjust(3))
            with self.subTest(name=name), self.assertRaises(ca.ComplexAnalysisError):
                ca.inspect_pdb(data)

    def test_multimodel_duplicate_and_alt(self):
        cases = [b'MODEL        1\nMODEL        2\n' + pdb(), pdb() + pdb(),
                 (atom(9, 'CA', 'ALA', 1, 'C', 'ATOM', alt='A') + '\n').encode() + pdb(),
                 (atom(9, 'CA', 'ALA', 1, 'C', 'ATOM', insertion='A') + '\n').encode() + pdb()]
        for data in cases:
            with self.subTest(data=data[:30]), self.assertRaises(ca.ComplexAnalysisError):
                ca.inspect_pdb(data)

    def test_bad_coordinates_and_input(self):
        for data in [b'', b'not pdb', b'\xff', b'X' * (ca.MAX_BYTES + 1), pdb().replace(b'   1.200', b'     nan')]:
            with self.subTest(length=len(data)), self.assertRaises(ca.ComplexAnalysisError):
                ca.inspect_pdb(data)


class ReportTests(unittest.TestCase):
    def test_hydrogen_donor_direction(self):
        result = ca.parse_plip_report(xml(), 'LIG:A:2')
        self.assertEqual(result['contacts'][0]['ligand_atoms'], [4])
        reversed_report = xml(extra='<protisdon>False</protisdon><acceptoridx>1</acceptoridx><donoridx>4</donoridx>')
        self.assertEqual(ca.parse_plip_report(reversed_report, 'LIG:A:2')['contacts'][0]['ligand_atoms'], [4])
        self.assertEqual(result['mapping'], {2: 0, 3: 1, 4: 2})

    def test_all_interaction_types(self):
        for kind, field, extra in [
            ('hydrophobic_interaction', 'dist', '<ligcarbonidx>2</ligcarbonidx>'),
            ('salt_bridge', 'dist', '<lig_idx_list><idx>3</idx><idx>4</idx></lig_idx_list>'),
            ('pi_stack', 'centdist', '<lig_idx_list><idx>2</idx></lig_idx_list>'),
            ('pi_cation_interaction', 'dist', '<lig_idx_list><idx>2</idx></lig_idx_list>')]:
            with self.subTest(kind=kind):
                result = ca.parse_plip_report(xml(kind, extra, field), 'LIG:A:2')
                self.assertEqual(result['contacts'][0]['kind'], kind)
                self.assertEqual(result['contacts'][0]['distance'], 3.1)

    def test_no_contacts_is_valid(self):
        result = ca.parse_plip_report(xml(with_contact=False), 'LIG:A:2')
        self.assertEqual(result['contacts'], [])
        self.assertIn('Расстояние', ca.contacts_csv(result).decode('utf-8-sig'))

    def test_invalid_reports(self):
        for report in ['bad xml', xml().replace('3.1', 'nan'), xml().replace('<composite>False', '<composite>True'),
                       xml().replace('<acceptoridx>4</acceptoridx>', '')]:
            with self.subTest(report=report[:25]), self.assertRaises(ca.ComplexAnalysisError):
                ca.parse_plip_report(report, 'LIG:A:2')
        with self.assertRaises(ca.ComplexAnalysisError):
            ca.parse_plip_report(xml(), 'OTHER:A:2')

    def test_csv_and_formula_escaping(self):
        result = ca.parse_plip_report(xml(), 'LIG:A:2')
        result['contacts'][0]['residue'] = '=unsafe'
        exported = ca.contacts_csv(result).decode('utf-8-sig')
        self.assertIn("'=unsafe", exported)
        self.assertIn('Донор — акцептор', exported)


class RunnerTests(unittest.TestCase):
    def test_isolated_run_and_cleanup(self):
        folders = []
        def run(command, **kwargs):
            folders.append(kwargs['cwd'])
            self.assertNotIn('-i', command)
            self.assertIn('--nofix', command)
            self.assertEqual(kwargs['timeout'], 90)
            Path(kwargs['cwd'], 'report.xml').write_text(xml(), encoding='utf-8')
            return types.SimpleNamespace(returncode=0, stderr='')
        with patch.object(ca.subprocess, 'run', side_effect=run):
            first = ca.analyze_complex(pdb(), 'LIG:A:2')
            ca.analyze_complex(pdb(), 'LIG:A:2')
        self.assertEqual(set(first['atoms']), {2, 3, 4})
        self.assertNotEqual(folders[0], folders[1])
        self.assertTrue(all(not Path(folder).exists() for folder in folders))

    def test_timeout_and_dependency_failure(self):
        with patch.object(ca.subprocess, 'run', side_effect=subprocess.TimeoutExpired('plip', 90)):
            with self.assertRaisesRegex(ca.ComplexAnalysisError, '90'):
                ca.analyze_complex(pdb(), 'LIG:A:2')
        with patch.object(ca.subprocess, 'run', return_value=types.SimpleNamespace(returncode=1, stderr='No module named plip')):
            with self.assertRaisesRegex(ca.ComplexAnalysisError, 'Open Babel'):
                ca.analyze_complex(pdb(), 'LIG:A:2')

    def test_wrong_ligand_never_runs(self):
        with patch.object(ca.subprocess, 'run') as run:
            with self.assertRaises(ca.ComplexAnalysisError):
                ca.analyze_complex(pdb(), 'LIG:B:99')
            run.assert_not_called()


class DiagramTests(unittest.TestCase):
    def test_unique_residue_layout(self):
        contact = dict(residue='ALA', chain='A', number=1)
        width, height, positions = diagram_layout([contact, contact, dict(contact, chain='B')])
        self.assertEqual(len(positions), 2)
        self.assertEqual(width, 1200)
        self.assertGreaterEqual(height, 740)

    @unittest.skipUnless(importlib.util.find_spec('rdkit'), 'RDKit is not installed in this test environment')
    def test_real_svg_png_and_mapping_rejection(self):
        from core.interaction_diagram import render_diagram
        result = ca.parse_plip_report(xml(), 'LIG:A:2')
        result['atoms'] = {k: a for k, a in ca.inspect_pdb(pdb())['atoms'].items() if k != 1}
        svg, png = render_diagram(result)
        self.assertIn(b'<svg', svg)
        self.assertTrue(png.startswith(b'\x89PNG'))
        result['mapping'] = {}
        with self.assertRaises(ca.ComplexAnalysisError):
            render_diagram(result)


if __name__ == '__main__':
    unittest.main()
