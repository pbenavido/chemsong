from util.reaction import molecules_to_bond_energy_df


def test_molecules_to_df():
    mols = {}
    mols[1] = ["CF", "CO"]
    mols[2] = ["c1ccccc1", "OC"]
    mols[3] = ["CCc1ccccc1"]

    df = molecules_to_bond_energy_df(mols)

    assert True
