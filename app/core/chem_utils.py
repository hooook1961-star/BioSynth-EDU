from rdkit import Chem
from rdkit.Chem import AllChem
import pubchempy as pcp

def smiles_to_3d_block(smiles: str, optimize=False) -> str:
    """Генерирует 3D-молекулу с опциональной минимизацией энергии."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        # Добавляем неявные водороды (важно для правильной геометрии)
        mol = Chem.AddHs(mol)
        
        # Генерируем начальные 3D координаты
        # ETKDG - современный стандарт генерации конформаций
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        if optimize:
            # Используем силовое поле MMFF94 для минимизации энергии
            AllChem.MMFFOptimizeMolecule(mol)
            
        return Chem.MolToMolBlock(mol)
    except Exception:
        return None

def get_pubchem_data(smiles: str):
    """Запрос физико-химических констант."""
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        if compounds:
            cmp = compounds[0]
            return {
                "mw": cmp.molecular_weight,
                "logp": cmp.xlogp,
                "formula": cmp.molecular_formula
            }
    except Exception:
        return None
