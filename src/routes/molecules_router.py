from fastapi import APIRouter, HTTPException, Depends, Query
from sqlalchemy.ext.asyncio import AsyncSession
from typing import List
from src import crud, schemas
from src.database import get_db
from src.logger import logger

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
    logger.info(f"Attempting to add molecule: {molecule}")
    try:
        result = await crud.create_molecule(db=db, molecule=molecule)
        logger.info(f"Molecule added successfully: {result}")
        return result
    except ValueError as ve:
        logger.error(f"Error adding molecule: {ve}")
        raise HTTPException(status_code=400, detail=str(ve))


@router.get("/", response_model=List[schemas.Molecule])
async def list_molecules(
    skip: int = Query(0, alias="start", ge=0),
    limit: int = Query(100, le=1000),  # You can adjust the max limit as needed
    db: AsyncSession = Depends(get_db)
):
    """
    List molecules from the database with pagination.

    :param skip: Number of records to skip (for pagination).
    :param limit: Maximum number of records to retrieve.
    :param db: The database session.
    :return: A list of molecules.
    """
    logger.info("Listing molecules with skip=%d and limit=%d", skip, limit)

    # Initialize an empty list to collect molecules
    molecules = []

    # Collect results from the asynchronous generator
    async for molecule in crud.get_molecules(db=db, skip=skip, limit=limit):
        molecules.append(molecule)

    return molecules


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
    logger.info(f"Fetching molecule with ID: {molecule_id}")
    molecule = await crud.get_molecule(db=db, molecule_id=molecule_id)
    if not molecule:
        logger.warning(f"Molecule with ID {molecule_id} not found.")
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
    logger.info(f"Updating molecule with ID: {molecule_id}")
    try:
        result = await crud.update_molecule(db=db, molecule_id=molecule_id, molecule=updated_molecule)
        logger.info(f"Updating molecule with ID: {molecule_id}")
        return result
    except ValueError as ve:
        logger.error(f"Error updating molecule with ID {molecule_id}: {ve}")
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
    logger.info(f"Deleting molecule with ID: {molecule_id}")
    try:
        result = await crud.delete_molecule(db=db, molecule_id=molecule_id)
        logger.info(f"Molecule with ID {molecule_id} deleted successfully.")
        return result
    except ValueError:
        logger.warning(f"Molecule with ID {molecule_id} not found for deletion.")
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
    logger.info(f"Searching for molecules containing substructure: {mol}")
    try:
        if not mol:
            logger.info(f"No molecules found containing substructure: {mol}")
            raise HTTPException(status_code=400,
                                detail="Invalid substructure SMILES string. Must be a non-empty string.")

        matches = await crud.substructure_search(db=db, mol=mol)

        if not matches:
            raise HTTPException(status_code=404, detail="No molecules found containing the given substructure.")
        logger.info(f"Found molecules containing substructure: {mol}")
        return {"matches": matches}
    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        logger.error(f"Value error during substructure search: {e}")
        raise HTTPException(status_code=500, detail=f"An unexpected error occurred: {str(e)}")
