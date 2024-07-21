from rdkit import Chem
from fastapi import FastAPI, HTTPException, File, UploadFile
from pydantic import BaseModel
from rdkit.Chem import Draw
import os
from fastapi.responses import FileResponse


app = FastAPI()
# Creating a dir for image storing, unrelated to base code, just for optional file uploading task.
IMAGE_DIR = "imgs"
os.makedirs(IMAGE_DIR, exist_ok=True)

smiles_db = [
    {"id": 1, "name": "Water", "smiles": "O"},
    {"id": 2, "name": "Benzene", "smiles": "c1ccccc1"},
    {"id": 3, "name": "Acetic acid", "smiles": "CC(=O)O"},
    {"id": 4, "name": "Aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
    {"id": 5, "name": "Ethanol", "smiles": "CCO"},
    {"id": 6, "name": "Methanol", "smiles": "CO"},
    {"id": 7, "name": "Propane", "smiles": "CCC"},
    {"id": 8, "name": "Methane", "smiles": "C"},
    {"id": 9, "name": "Ethane", "smiles": "CC"}
]


class Molecule(BaseModel):
    id: int
    name: str
    smiles: str


@app.post("/molecules/add", response_model=Molecule, response_description="Added molecule", status_code=201)
def add_molecule(molecule: Molecule):
    """
    **Add a new molecule to the database.**

    **Parameters**:
      - **molecule** (*according to Class "Molecule"-s structure*): The molecule to add.

    **Returns**:
      - **Molecule**: The added molecule.

    **Raises**:
      - **HTTPException**: If a molecule with the same ID or name already exists.
    """
    for mol in smiles_db:
        if mol["id"] == molecule.id:
            raise HTTPException(status_code=400, detail="Molecule with this ID already exists.")
        if mol["name"].lower() == molecule.name.lower():
            raise HTTPException(status_code=400, detail="Molecule with this name already exists.")
    # Using model_dump() instead of dict() because Pydantic deprecated it.
    smiles_db.append(molecule.model_dump())
    return molecule


@app.get("/molecules/", response_description="List of all molecules")
def list_molecules():
    """
      **List all molecules.**

      **Returns**:
        - **List[Molecule]**: A list of all molecules.
      """
    return smiles_db


@app.get("/molecules/{molecule_id}", response_description="Returned Molecule by ID")
def get_molecule_by_id(molecule_id: int):
    """
    **Get a molecule by its ID.**

    **Parameters**:
      - **molecule_id** (*int*): The ID of the molecule to retrieve.

    **Returns**:
      - **Molecule**: The molecule with the given ID.

    **Raises**:
      - **HTTPException**: If the molecule is not found.
    """
    for mol in smiles_db:
        if mol["id"] == molecule_id:
            return mol
    raise HTTPException(status_code=404, detail="Molecule not found.")


@app.put("/molecules/{molecule_id}", response_model=Molecule, response_description="Updated molecule")
def update_molecule_by_id(molecule_id: int, updated_molecule: Molecule):
    """
       **Update an existing molecule by its ID.**

       **Parameters**:
         - **molecule_id** (*int*): The ID of the molecule to update.
         - **updated_molecule** (*Molecule*): The updated molecule data.

       **Returns**:
         - **Molecule**: The updated molecule.

       **Raises**:
         - **HTTPException**: If the molecule is not found.
       """
    for index, mol in enumerate(smiles_db):
        if mol["id"] == molecule_id:
            smiles_db[index] = updated_molecule.model_dump()
            return updated_molecule.model_dump()
    raise HTTPException(status_code=404, detail="Molecule not found.")


@app.delete("/molecules/{molecule_id}", response_description="Deleted molecule",)
def delete_molecule_by_id(molecule_id: int):
    """
       **Delete a molecule by its ID.**

       **Parameters**:
         - **molecule_id** (*int*): The ID of the molecule to delete.

       **Returns**:
         - **Deleted Molecule**: Returns the removed molecule and its data.

       **Raises**:
         - **HTTPException**: If the molecule is not found.
       """
    for index, mol in enumerate(smiles_db):
        if mol["id"] == molecule_id:
            deleted_user = smiles_db.pop(index)
            return deleted_user
    raise HTTPException(status_code=404, detail="Molecule not found.")


@app.get("/molecules/search/", response_description="List of all molecules containing the substructure")
def search_substructure(mol: str):
    """
    **Search for molecules containing the given substructure.**

    **Parameters**:
      - **mol** (*str*): SMILES string representing the substructure to search for.

    **Returns**:
      - **Dict[str, List[str]]**: List of SMILES strings of molecules containing the substructure.
    """
    try:
        matches = substructure_search([mol_dict["smiles"] for mol_dict in smiles_db], mol)
    except ValueError as ve:
        # Provided error's desc from base function
        raise HTTPException(status_code=400, detail=str(ve))
    return {"matches": matches}


def substructure_search(mols, mol):
    """
      Func returns a list of all molecules from argument 1 that contain substructure from argument 2.

      Arguments:
          mols (list): A list of SMILES strings containing the molecules to search.
          mol (str): A SMILES string which has the substructure to search for.

      Returns:
          list: A list of SMILES strings representing the molecules that contain the substructure.
          str: An error message if an exception occurs.
      """
    substructure = Chem.MolFromSmiles(mol)
    if substructure is None:
        raise ValueError("Unknown substructure SMILES string")

    matches = []
    for mol_smiles in mols:
        try:
            molecule = Chem.MolFromSmiles(mol_smiles)
            if molecule is None:
                raise ValueError(f"Invalid molecule SMILES string: {mol_smiles}")

            if molecule.HasSubstructMatch(substructure):
                matches.append(mol_smiles)
        except ValueError as ve:
            print(ve)

    return matches


# Checking substructure function and it's correctness
final = substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
expected_result = ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
assert final == expected_result, f"Expected {expected_result}, got {final}"


# Experimenting ->
@app.post("/molecules/upload/", response_description="Succesfully created images for the molecules")
async def upload_file(file: UploadFile = File(...)):
    """
    **Upload a text file containing SMILES strings and generate images for each molecule.**

    **Parameters**:
      - **file** (*UploadFile*): A text file (`.txt`) containing SMILES strings, each on a new line.

    **Returns**:
      - **JSON**:  JSON  with a list of molecules, where each entry includes:
        - **image_url** (*str*): The URL to the generated image of the molecule.
        - **smiles** (*str*): The SMILES string used to generate the molecule image.

    **Raises**:
      - **HTTPException**:
        - **400 Bad Request**: If the file type is not `.txt` or if any SMILES string is invalid.
    """
    if not file.filename.endswith(".txt"):
        raise HTTPException(status_code=400, detail="Invalid file type. Please upload a .txt file.")

    contents = await file.read()
    smiles_list = contents.decode().splitlines()

    results = []

    for i, smiles in enumerate(smiles_list):
        try:
            image_path = os.path.join(IMAGE_DIR, f"molecule_{i}.png")
            smiles = smiles.strip()
            if smiles:  # Making sure that it is not an empty string
                draw_molecule(smiles, image_path)
                results.append({
                    "image_url": image_path,
                    "smiles": smiles
                })
        except ValueError as e:
            return HTTPException(status_code=400, detail=str(e))

    return {"molecules": results}


def draw_molecule(smiles: str, image_path: str):
    """
    Draw a molecule from a SMILES string and save it as an image.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    img = Draw.MolToImage(mol)
    img.save(image_path)


@app.get("/molecules/images/{image_name}", response_description="Image of the molecule")
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
    image_path = os.path.join(IMAGE_DIR, image_name)
    if not os.path.isfile(image_path):
        raise HTTPException(status_code=404, detail="Image not found")
    return FileResponse(image_path)
