from rdkit import Chem
from rdkit.Chem import Draw

from loguru import logger

def render(mols):
    for step, smiles_list in mols.items():
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
        img = Draw.MolsToGridImage(mol_list, molsPerRow=2, subImgSize=(200, 200))
        img.show()
        logger.info(f"Displayed molecules for step {step}")