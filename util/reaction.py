import pandas as pd
from loguru import logger

from util.bond_energies import bond_energies, list_bonds


def molecules_to_bond_energy_df(mols):

    step_bond_energies = {}  # Dictionary to store bond energies for each step

    for stepnum in mols:
        logger.info(f"step {stepnum}: {mols[stepnum]}")
        bond_energies_list = []
        for mol in mols[stepnum]:
            bonds = list_bonds(mol)
            logger.info(f"The compound {mol} contains the following bonds: {bonds}")

            for bond in bonds:
                bond_energy = bond_energies.get(bond, None)
                if bond_energy is not None:
                    bond_energies_list.append(bond_energy)
                else:
                    logger.warning(f"Energy for bond {bond} not found. Skipping.")

        step_bond_energies[stepnum] = bond_energies_list

    # Creating a DataFrame from the dictionary
    bond_df = pd.DataFrame(
        list(step_bond_energies.items()), columns=["Step", "Bond Energies"]
    )
    logger.success(f"\n{bond_df}")

    return bond_df
