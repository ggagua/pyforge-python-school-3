import io
from rdkit.Chem import Draw


# In utils since has reusable value
def draw_molecule(mol) -> bytes:
    """
    Draw a molecule from a SMILES string and return the image as bytes.

    Parameters:
        mol (Mol): The RDKit molecule object to draw.

    Returns:
        bytes: The image of the molecule in PNG format.
    """
    img = Draw.MolToImage(mol)
    buf = io.BytesIO()
    img.save(buf, format='PNG')
    buf.seek(0)
    return buf.getvalue()
