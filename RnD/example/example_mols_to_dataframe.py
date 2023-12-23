import pandas as pd

# Typical bond energies (kJ/mol, approximate values from literature)
bond_energies = {
    "Bond Type": ["C-C (aromatic)", "C-C (aliphatic)", "C-H", "C-Cl"],
    "Energy (kJ/mol)": [615, 348, 414, 328],  # Example values
}

df = pd.DataFrame(bond_energies)
print(df)
