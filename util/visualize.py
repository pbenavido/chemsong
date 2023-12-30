import io

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw


def verify_smiles(smiles):
    return Chem.MolFromSmiles(smiles) is not None


def render_smiles(smiles_list):
    mol_list = [
        Chem.MolFromSmiles(smile) for smile in smiles_list if verify_smiles(smile)
    ]
    # Convert RDKit image to PIL Image
    img = Draw.MolsToGridImage(
        mol_list, molsPerRow=4, subImgSize=(300, 100), useSVG=False
    )
    img_byte_arr = io.BytesIO()
    img.save(img_byte_arr, format="PNG")
    img_byte_arr = img_byte_arr.getvalue()
    pil_img = Image.open(io.BytesIO(img_byte_arr))
    return pil_img
