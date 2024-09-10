# import pytest
# from httpx import AsyncClient, ASGITransport
# from sqlalchemy.ext.asyncio import AsyncSession
# from src.models import Molecule
# from src.schemas import MoleculeCreate
# from src.main import app
#
# @pytest.mark.asyncio
# async def test_add_molecule():
#     async with AsyncClient(
#             transport=ASGITransport(app=app), base_url="http://test"
#     ) as ac:
#         response = await ac.post(
#             "/molecules/add",
#             json={
#                 "name": "Test Moleculefd",
#                 "smiles": "C1=CC=CC=C1dfa",  # Example SMILES string
#             }
#         #
#         assert response.status_code == 200
#         assert response.json()["name"] == "Test Moleculefd"
#         assert response.json()["smiles"] == "C1=CC=CC=C1dfa"
