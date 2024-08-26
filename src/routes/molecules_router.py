from fastapi import APIRouter, HTTPException, Depends
from sqlalchemy.ext.asyncio import AsyncSession
from typing import List
from src import crud, schemas
from src.database import get_db

router = APIRouter()


@router.post("/add", response_model=schemas.Molecule)
async def add_molecule(molecule: schemas.MoleculeCreate, db: AsyncSession = Depends(get_db)):
    """
    Add a new molecule to the database.

    **Parameters**:
      - **molecule**: The data required to create a new molecule.

    **Returns**:
      - **Molecule**: The created molecule with its details.

    **Raises**:
      - **HTTPException**: If a molecule with the same ID or name already exists.
    """
    try:
        return await crud.create_molecule(db=db, molecule=molecule)
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))


@router.get("/", response_model=List[schemas.Molecule])
async def list_molecules(db: AsyncSession = Depends(get_db)):
    """
    List all molecules in the database.

    **Returns**:
      - **List[schemas.Molecule]**: A list of all molecules with their details.
    """
    return await crud.get_molecules(db=db)


@router.get("/{molecule_id}", response_model=schemas.Molecule)
async def get_molecule_by_id(molecule_id: int, db: AsyncSession = Depends(get_db)):
    """
    Get details of a specific molecule by its ID.

    **Parameters**:
      - **molecule_id** (*int*): The ID of the molecule to retrieve.

    **Returns**:
      - **Molecule**: The molecule with the specified ID.

    **Raises**:
      - **HTTPException**: If the molecule is not found.
    """
    molecule = await crud.get_molecule(db=db, molecule_id=molecule_id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found.")
    return molecule


@router.put("/{molecule_id}", response_model=schemas.Molecule)
async def update_molecule_by_id(molecule_id: int, updated_molecule: schemas.MoleculeUpdate,
                                db: AsyncSession = Depends(get_db)):
    """
    Update details of a specific molecule by its ID.

    **Parameters**:
      - **molecule_id** (*int*): The ID of the molecule to update.
      - **updated_molecule** (*schemas.MoleculeUpdate*): The updated data for the molecule.

    **Returns**:
      - **Molecule**: The updated molecule.

    **Raises**:
      - **HTTPException**: If the molecule is not found.
    """
    try:
        return await crud.update_molecule(db=db, molecule_id=molecule_id, molecule=updated_molecule)
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))


@router.delete("/{molecule_id}", response_model=schemas.Molecule)
async def delete_molecule_by_id(molecule_id: int, db: AsyncSession = Depends(get_db)):
    """
    Delete a specific molecule by its ID.

    **Parameters**:
      - **molecule_id** (*int*): The ID of the molecule to delete.

    **Returns**:
      - **Molecule**: The deleted molecule.

    **Raises**:
      - **HTTPException**: If the molecule is not found.
    """
    try:
        return await crud.delete_molecule(db=db, molecule_id=molecule_id)
    except ValueError:
        raise HTTPException(status_code=404, detail="Molecule not found.")


@router.get("/molecules/search/", response_description="Returned the list of all molecules containing the substructure",
            summary="Search for molecules containing a substructure", tags=["Molecules"])
async def search_molecules_by_substructure(mol: str, db: AsyncSession = Depends(get_db)) -> dict:
    """
    Search for molecules containing the given substructure.

    **Parameters**:
      - **mol** (*str*): SMILES string representing the substructure to search for.

    **Returns**:
      - **dict**: A dictionary containing a list of SMILES strings of molecules containing the substructure.

    **Raises**:
      - **HTTPException**: If the substructure SMILES string is invalid or if there are issues processing the molecules.
    """
    try:
        if not mol:
            raise HTTPException(status_code=400,
                                detail="Invalid substructure SMILES string. Must be a non-empty string.")

        matches = await crud.substructure_search(db=db, mol=mol)

        if not matches:
            raise HTTPException(status_code=404, detail="No molecules found containing the given substructure.")

        return {"matches": matches}
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"An unexpected error occurred: {str(e)}")
