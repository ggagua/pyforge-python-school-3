from celery_worker import celery
from crud import substructure_search
from database import get_db
import asyncio


@celery.task
def substructure_search_task(mol: str):
    """
    Perform substructure search as a Celery task.

    :param mol: The SMILES string representing the substructure.
    :return: A list of SMILES strings for molecules containing the substructure.
    """

    async def async_search():
        async for db in get_db():  #
            return await substructure_search(db, mol)

    loop = asyncio.get_event_loop()
    return loop.run_until_complete(async_search())