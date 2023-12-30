from loguru import logger
from rdkit import Chem

# bond energy values in Kj/mol at 273K
# transcribed from source: Data from J. E. Huheey, E. A. Keiter, and R. L. Keiter, Inorganic Chemistry, 4th ed. (1993).
bond_energies = {
    "O-O": 142,
    "O-O": 142,
    "I-I": 149,
    "I-I": 149,
    "F-F": 155,
    "F-F": 155,
    "N-N": 167,
    "N-N": 167,
    "Br-I": 175,
    "I-Br": 175,
    "O-F": 190,
    "F-O": 190,
    "Br-Br": 190,
    "Br-Br": 190,
    "N-O": 201,
    "O-N": 201,
    "O-Br": 201,
    "Br-O": 201,
    "O-I": 201,
    "I-O": 201,
    "P-P": 201,
    "P-P": 201,
    "Cl-I": 208,
    "I-Cl": 208,
    "C-I": 213,
    "I-C": 213,
    "Cl-Br": 216,
    "Br-Cl": 216,
    "O-Cl": 218,
    "Cl-O": 218,
    "S-Br": 218,
    "Br-S": 218,
    "Si-Si": 222,
    "Si-Si": 222,
    "S-S": 226,
    "S-S": 226,
    "Cl-Cl": 240,
    "Cl-Cl": 240,
    "N-Br": 243,
    "Br-N": 243,
    "F-Cl": 249,
    "Cl-F": 249,
    "F-Br": 249,
    "Br-F": 249,
    "S-Cl": 255,
    "Cl-S": 255,
    "C-S": 272,
    "S-C": 272,
    "F-I": 278,
    "I-F": 278,
    "N-F": 283,
    "F-N": 283,
    "S-F": 284,
    "F-S": 284,
    "C-Br": 285,
    "Br-C": 285,
    "H-I": 295,
    "I-H": 295,
    "C-N": 305,
    "N-C": 305,
    "N-Cl": 313,
    "Cl-N": 313,
    "C-Si": 318,
    "Si-C": 318,
    "H-Si": 318,
    "Si-H": 318,
    "H-P": 322,
    "P-H": 322,
    "C-Cl": 327,
    "Cl-C": 327,
    "C-C": 346,
    "C-C": 346,
    "C-O": 358,
    "O-C": 358,
    "H-Br": 362,
    "Br-H": 362,
    "H-S": 363,
    "S-H": 363,
    "H-N": 386,
    "N-H": 386,
    "H-C": 411,
    "C-H": 411,
    "N=N": 418,
    "N=N": 418,
    "H-Cl": 428,
    "Cl-H": 428,
    "H-H": 432,
    "H-H": 432,
    "Si-O": 452,
    "O-Si": 452,
    "H-O": 459,
    "O-H": 459,
    "C-F": 485,
    "F-C": 485,
    "O=O": 494,
    "O=O": 494,
    "S=O": 532,
    "O=S": 532,
    "H-F": 565,
    "F-H": 565,
    "C=C": 602,
    "C=C": 602,
    "N=O": 607,
    "O=N": 607,
    "C=N": 615,
    "N=C": 615,
    "C=O": 749,
    "O=C": 749,
    "C=:C": 835,
    "C=:C": 835,
    "C=:N": 887,
    "N=:C": 887,
    "N=:N": 942,
    "N=:N": 942,
    "C=:O": 1072,
    "O=:C": 1072,
}
# there are absolutely better ways to encode that dictionary, I'm so sorry


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
