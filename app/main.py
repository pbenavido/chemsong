'''
chemsong main
run with: python -m app.main
'''
from util.reaction import molecules_to_bond_energy_df
from util.visualize import *
from util.midi import df_to_notes



def chemsong():
    # TODO: get mols from user input
    mols = {}
    mols[1] = ["c1ccccc1"]
    mols[2] = ["c1ccccc1", "CC"]
    mols[3] = ["CCc1ccccc1"]

    # get bond energies from mols
    bond_df = molecules_to_bond_energy_df(mols)

    # visualize steps
    render(mols)

    # play notes
    df_to_notes(bond_df)

    # TODO: sync notes and visualizations

chemsong()