from rdkit import Chem
from loguru import logger

# bond energy values in Kj/mol at 273K
# transcribed from source: Data from J. E. Huheey, E. A. Keiter, and R. L. Keiter, Inorganic Chemistry, 4th ed. (1993).
bond_energies = {
"O-O": 142,
"I-I": 149,
"F-F": 155,
"N-N": 167,
"Br-I": 175,
"O-F": 190,
"Br-Br": 190,
"N-O": 201,
"O-Br": 201,
"O-I": 201,
"P-P": 201,
"Cl-I": 208,
"C-I": 213,
"Cl-Br": 216,
"O-Cl": 218,
"S-Br": 218,
"Si-Si": 222,
"S-S": 226,
"Cl-Cl": 240,
"N-Br": 243,
"F-Cl": 249,
"F-Br": 249,
"S-Cl": 255,
"C-S": 272,
"F-I": 278,
"N-F": 283,
"S-F": 284,
"C-Br": 285,
"H-I": 295,
"C-N": 305,
"N-Cl": 313,
"C-Si": 318,
"H-Si": 318,
"H-P": 322,
"C-Cl": 327,
"C-C": 346,
"C-O": 358,
"H-Br": 362,
"H-S": 363,
"H-N": 386,
"H-C": 411,
"N=N": 418,
"H-Cl": 428,
"H-H": 432,
"Si-O": 452,
"H-O": 459,
"C-F": 485,
"O=O": 494,
"S=O": 532,
"H-F": 565,
"C=C": 602,
"N=O": 607,
"C=N": 615,
"C=O": 749,
"C=:C": 835,
"C=:N": 887,
"N=:N": 942,
"C=:O": 1072,
}

def list_bonds(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        bonds = []
        for bond in mol.GetBonds():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            bond_type = bond.GetBondTypeAsDouble()
            bond_symbol = f"{begin_atom.GetSymbol()}-{end_atom.GetSymbol()}"
            if bond_type == 1.5:  # Handling aromatic bonds
                bond_symbol = f"{begin_atom.GetSymbol()}=:{end_atom.GetSymbol()}"
            elif bond_type == 2.0:
                bond_symbol = f"{begin_atom.GetSymbol()}={end_atom.GetSymbol()}"
            elif bond_type == 3.0:
                bond_symbol = f"{begin_atom.GetSymbol()}=:{end_atom.GetSymbol()}"
            bonds.append(bond_symbol)
        return list(set(bonds))  # Remove duplicates
    except Exception as e:
        logger.error(f"Invalid SMILES string {smiles}: {e}")
        return []
