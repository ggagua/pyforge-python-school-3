from fastapi import APIRouter, HTTPException, UploadFile, File
from fastapi.responses import StreamingResponse
from src.utils import draw_molecule
from rdkit import Chem
import io

# Global storage to store images temporarily
image_storage = {}

router = APIRouter()


@router.post("/upload/", response_description="Successfully created images for the molecules", status_code=201,
             summary="Upload a text file containing SMILES strings, parse the file and generate images for each molecule",
             tags=["Files"])
async def upload_file(file: UploadFile = File(...)):
    """
    **Upload a text file containing SMILES strings and generate images for each molecule**
    ***(only valid SMILES will be processed and rest will be ignored).***

    **Parameters**:
      - **file** (*UploadFile*): A text file (`.txt`) containing SMILES strings, each on a new line.

    **Returns**:
      - **JSON**: JSON with a list of molecules, where each entry includes:
        - **image_name** (*str*): The name of the generated image of the molecule.
        - **smiles** (*str*): The SMILES string used to generate the molecule image.

    **Raises**:
      - **HTTPException**:
        - **400 Bad Request**: If the file type is not `.txt`.
        - **422 Unprocessable Entity**: If every single data in the file is invalid (not SMILES).
    """
    if not file.filename.endswith(".txt"):
        raise HTTPException(status_code=400, detail="Invalid file type. Please upload a .txt file.")

    contents = await file.read()
    smiles_list = contents.decode().split()

    results = []

    for i, smiles in enumerate(smiles_list):
        try:
            smiles = smiles.strip()
            if smiles:  # Making sure that it is not an empty string
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:  # If the given SMILES is valid, include it in output
                    image_name = f"molecule_{i}.png"
                    image_data = draw_molecule(mol)
                    image_storage[image_name] = image_data  # Store the image with its name
                    results.append({
                        "image_name": image_name,
                        "smiles": smiles
                    })  # Return name and SMILES instead of returning the original dict
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))
    if not results:
        raise HTTPException(status_code=422, detail="File did not contain a single valid SMILES string.")
    return {"molecules": results}


@router.get("/images/{image_name}", response_description="Image of the molecule",
            summary="Get an image by its name",
            tags=["Files"])
async def get_image(image_name: str):
    """
    **Return an image of a molecule by its image name.**

    **Parameters**:
      - **image_name** (*str*): The name of the image file to get.

    **Returns**:
      - **File**: The image file of the requested molecule.

    **Raises**:
      - **HTTPException**: If the image file is not found.
    """
    image_data = image_storage.get(image_name)
    if not image_data:
        raise HTTPException(status_code=404, detail="Image not found")

    return StreamingResponse(io.BytesIO(image_data), media_type="image/png")
