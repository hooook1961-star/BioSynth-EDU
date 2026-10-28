"""Deployment smoke check: imports, chemistry bindings, drawing and PLIP CLI.

Run in the same Python environment as Streamlit after installing requirements.txt.
No network calls and no user PDB files are needed.
"""
from importlib.metadata import version
import subprocess
import sys


def main():
    try:
        from openbabel import pybel
        from plip.structure.preparation import PDBComplex
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        from PIL import Image
        import io

        assert pybel.readstring('smi', 'CCO').molwt > 0
        assert PDBComplex() is not None
        molecule = Chem.MolFromSmiles('CCO')
        drawer = rdMolDraw2D.MolDraw2DCairo(240, 180)
        drawer.DrawMolecule(molecule)
        drawer.FinishDrawing()
        Image.open(io.BytesIO(drawer.GetDrawingText())).verify()
        completed = subprocess.run([sys.executable, '-m', 'plip.plipcmd', '--help'],
                                   capture_output=True, text=True, timeout=30)
        if completed.returncode:
            raise RuntimeError(completed.stderr.strip())
        for package in ['plip', 'openbabel', 'rdkit', 'Pillow']:
            print(f'{package}: {version(package)}')
        print('Analysis dependencies, PLIP CLI and PNG drawing: OK')
        return 0
    except Exception as exc:
        print(f'Analysis runtime check failed: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
