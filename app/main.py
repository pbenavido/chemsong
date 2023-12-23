'''
chemsong main
run with: python -m app.main
'''

from loguru import logger

import pygame
import time

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw

from util.bond_energies import bond_energies
from util.reaction import molecules_to_bond_energy_df
from util.visualize import *
from util.midi import *

mols = {}
mols[1] = ["c1ccccc1"]
mols[2] = ["c1ccccc1", "CC"]
mols[3] = ["CCc1ccccc1"]

bond_df = molecules_to_bond_energy_df(mols)

# Assuming bond_energies is already sorted
sorted_bond_energies = sorted(bond_energies.items(), key=lambda item: item[1])
sorted_bond_energies_dict = {k: i for i, (k, _) in enumerate(sorted_bond_energies)}

# normalize energies for midi mapping
bond_df['Normalized Energies'] = bond_df['Bond Energies'].apply(
    lambda energies: [sorted_bond_energies_dict.get(e) for e in energies]
)

# begin midi code section
def get_note_frequency(position):
    # Define the Arabic minor scale frequencies (in Hz) for simplicity
    arabic_minor_scale = [261.63, 293.66, 311.13, 349.23, 392.00, 415.30, 466.16]  # C4 scale
    return arabic_minor_scale[position % len(arabic_minor_scale)]

# Initialize pygame mixer
pygame.mixer.init()

# Iterate over DataFrame and play notes
for _, row in bond_df.iterrows():
    for pos in row['Normalized Energies']:
        frequency = get_note_frequency(pos)
        # Play the note (assuming a simple square wave for demonstration)
        # You might need a more sophisticated approach for realistic sound synthesis
        pygame.mixer.Sound(frequency=frequency, size=-16, array=None, buffer=4096).play()
        time.sleep(0.5)  # Delay between notes
# end midi code section

# begin visualization code
for step, smiles_list in mols.items():
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    img = Draw.MolsToGridImage(mol_list, molsPerRow=2, subImgSize=(200, 200))
    img.show()
    logger.info(f"Displayed molecules for step {step}")
# end visualization code