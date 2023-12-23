from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# Create molecules
benzene = Chem.MolFromSmiles('c1ccccc1')
ethyl_chloride = Chem.MolFromSmiles('CCCl')  # Assuming ethyl chloride as the ethyl group donor
ethylbenzene = Chem.MolFromSmiles('CCc1ccccc1')

# Generate 2D coordinates for better visualization
AllChem.Compute2DCoords(benzene)
AllChem.Compute2DCoords(ethyl_chloride)
AllChem.Compute2DCoords(ethylbenzene)

# Visualize the molecules
img = Draw.MolsToGridImage([benzene, ethyl_chloride, ethylbenzene], 
                           molsPerRow=3, subImgSize=(200, 200), 
                           legends=['Benzene', 'Ethyl Chloride', 'Ethylbenzene'])
img.show()
