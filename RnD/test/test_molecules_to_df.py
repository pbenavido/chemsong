"""
WIP
"""

from loguru import logger

from util.reaction import molecules_to_bond_energy_df

mols = {}
mols[1] = ["c1ccccc1"]
mols[2] = ["c1ccccc1", "CC"]
mols[3] = ["CCc1ccccc1"]

df = molecules_to_bond_energy_df(mols)

breakpoint()

assert df
