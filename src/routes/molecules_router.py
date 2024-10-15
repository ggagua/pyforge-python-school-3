from fastapi import APIRouter, HTTPException, Depends, Query
from sqlalchemy.ext.asyncio import AsyncSession
from celery.result import AsyncResult
from src.tasks import substructure_search_task
from src.celery_worker import celery
from typing import List, Dict, Optional
from src import crud
from src.database import get_db
from src.logger import logger
import json
import redis
from src.schemas import Molecule, MoleculeCreate, MoleculeUpdate

router = APIRouter()

# If you are running this with docker, set the host to 'redis', if not - 'localhost'
redis_client = redis.Redis(host='redis', port=6379, db=0, decode_responses=True)


def get_cached_result(key: str) -> Optional[Dict]:
    """Retrieve cached result from Redis."""
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: Optional[Dict], expiration: int = 60):
    """Store result in Redis with expiration time."""
    if value is not None:
        redis_client.setex(key, expiration, json.dumps(value))


@router.post("/add", response_model=Molecule)
async def add_molecule(molecule: MoleculeCreate, db: AsyncSession = Depends(get_db)):
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


@router.get("/", response_model=List[Molecule])
async def list_molecules(
        skip: int = Query(0, alias="start", ge=0),
        limit: int = Query(100, le=1000),  # You can adjust the max limit as needed
        db: AsyncSession = Depends(get_db)
):
    """
    Get details of a specific molecule by its ID.
    **Parameters**:
      - **start** (*int*): Number from which you would like to fetch/start.
        - **limit** (*int*): Number of records to fetch.
    **Returns**:
      - **Molecule**: The molecule with the specified ID.
    **Raises**:
      - **HTTPException**: If the molecule is not found.
    """
    logger.info("Listing molecules with skip=%d and limit=%d", skip, limit)

    # Initialize an empty list to collect molecules
    molecules = []

    # Collect results from the asynchronous generator
    async for molecule in crud.get_molecules(db=db, skip=skip, limit=limit):
        molecules.append(molecule)

    return molecules


@router.get("/{molecule_id}", response_model=Molecule)
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


@router.put("/{molecule_id}", response_model=Molecule)
async def update_molecule_by_id(molecule_id: int, updated_molecule: MoleculeUpdate,
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


@router.delete("/{molecule_id}", response_model=Molecule)
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


@router.post("/tasks/search_molecules/")
async def start_substructure_search(mol: str):
    """
    Get the result of the substructure search task.

    **Parameters**:
      - **task_id** (*str*): The ID of the Celery task.

    **Returns**:
      - **object**: containing the task status and result, with the following keys:
        - **task_id** (*str*): The task ID.
        - **status** (*str*): The status of the task. It can be "Task is still processing", "Task completed", or state.
        - **result** (*any*, optional): The result of the task, available if the task is completed successfully.

    **Raises**:
      - **HTTPException**: If the task is not found or another error occurs.
    """
    if not mol:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES string. Must be a non-empty string.")

    task = substructure_search_task.delay(mol)
    return {"task_id": task.id}


@router.get("/tasks/{task_id}")
async def get_substructure_search_result(task_id: str):
    """
    Get the result of the substructure search task.

    **Parameters**:
      - **task_id** (*str*): The ID of the Celery task to query.

    **Returns**:
      - **dict**: A dictionary containing the status and result of the task. The dictionary includes:
        - **task_id** (*str*): The ID of the task being queried.
        - **status** (*str*): The current state of the task. Possible values are:
          - "Task is still processing" if the task is not yet completed.
          - "Task completed" if the task has finished successfully.
          - The task’s current state if it is neither pending nor successful.
        - **result** (*any*, optional): The result of the task if it is completed successfully.
    """
    task_result = AsyncResult(task_id, app=celery)

    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {"task_id": task_id, "status": "Task completed", "result": task_result.result}
    else:
        raise HTTPException(status_code=404, detail={"task_id": task_id, "status": task_result.state})
