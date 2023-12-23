from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# Create a benzene molecule
benzene = Chem.MolFromSmiles("c1ccccc1")

# Create a methyl-substituted benzene (toluene)
toluene = Chem.MolFromSmiles("Cc1ccccc1")

# Create a ethyl group
ethyl = Chem.MolFromSmiles("CC")

# Show ethyl group with benzene
"""
NOTE: Chem.CombineMols seems to just put the molecules next ot one another rather than combine the structures.
Need to specify bonding location/type?
"""
ethyl_and_benzene = Chem.CombineMols(benzene, ethyl)

# Create an ethyl-substituted benzene (ethylbenzene)
ethylbenzene = Chem.MolFromSmiles("CCc1ccccc1")

# Visualize the molecules
img = Draw.MolsToGridImage(
    [ethyl_and_benzene, ethylbenzene], molsPerRow=2, subImgSize=(300, 300)
)
img.save("experiment/rdkit_demo_image.png")

# # # Generate 3D coordinates
# AllChem.EmbedMolecule(hexane)
# AllChem.EmbedMolecule(methyl_hexane)
