## NOT WORKING! ##
## WILL BE FIXED LATER ##

# import pytest
# import httpx
# from sqlalchemy.ext.asyncio import AsyncSession
# from sqlalchemy.future import select
# from src.main import app
# from src.database import Base, get_db
# from src.models import Molecule
# from src.crud import substructure_search
#
# @pytest.fixture
# async def async_client():
#     async with httpx.AsyncClient(app=app, base_url="http://testserver") as client:
#         yield client
#
# @pytest.fixture
# async def db_session() -> AsyncSession:
#     async with get_db() as session:
#         async with session.begin():
#             await session.run_sync(Base.metadata.create_all)
#         try:
#             yield session
#         finally:
#             async with session.begin():
#                 await session.run_sync(Base.metadata.drop_all)
#
# @pytest.mark.parametrize(
#     "mols, substructure, expected",
#     [
#         (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1", ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),
#         (["CCO", "CCC", "CCN"], "C", ["CCO", "CCC", "CCN"]),
#         (["CCO", "CCC", "CCN"], "N", ["CCN"]),
#         (["CCO", "CCC", "CCN"], "O", ["CCO"]),
#         (["CCO", "CCC"], "CC", ["CCO", "CCC"]),
#     ]
# )
# @pytest.mark.asyncio
# async def test_substructure_search_param(mols, substructure, expected, db_session: AsyncSession):
#     # Clear the table
#     async with db_session() as session:
#         result = await session.execute(select(Molecule))
#         molecules = result.scalars().all()
#         for molecule in molecules:
#             await session.delete(molecule)
#         await session.commit()
#
#     # Add new molecules
#     async with db_session() as session:
#         for mol in mols:
#             session.add(Molecule(name=mol, smiles=mol))
#         await session.commit()
#
#     # Test the substructure search
#     async with db_session() as session:
#         result = await substructure_search(session, substructure)
#         assert result == expected
#
# @pytest.mark.asyncio
# async def test_substructure_search_invalid_substructure(db_session: AsyncSession):
#     # Clear the table
#     async with db_session() as session:
#         result = await session.execute(select(Molecule))
#         molecules = result.scalars().all()
#         for molecule in molecules:
#             await session.delete(molecule)
#         await session.commit()
#
#     # Test with an invalid substructure
#     async with db_session() as session:
#         with pytest.raises(ValueError, match="Unknown substructure SMILES string"):
#             await substructure_search(session, "invalid_smiles")
#
# @pytest.mark.asyncio
# async def test_substructure_search_invalid_molecule_smiles(db_session: AsyncSession):
#     # Clear the table
#     async with db_session() as session:
#         result = await session.execute(select(Molecule))
#         molecules = result.scalars().all()
#         for molecule in molecules:
#             await session.delete(molecule)
#         await session.commit()
#
#     # Add molecules including an invalid one
#     async with db_session() as session:
#         for mol in ["CCO", "CCC", "invalid_smiles"]:
#             session.add(Molecule(name=mol, smiles=mol))
#         await session.commit()
#
#     # Test with valid substructure but invalid molecule SMILES
#     async with db_session() as session:
#         with pytest.raises(ValueError, match="Invalid molecule SMILES string: invalid_smiles"):
#             await substructure_search(session, "CC")
#
# @pytest.mark.asyncio
# async def test_get_server(async_client):
#     response = await async_client.get("/")
#     assert response.status_code == 200
#     assert "server_id" in response.json()
#
#
# @pytest.mark.asyncio
# async def test_add_molecule(async_client):
#     new_molecule = {
#         "id": 10,
#         "name": "Caffeine",
#         "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
#     }
#     response = await async_client.post("/molecules/add", json=new_molecule)
#     assert response.status_code == 201
#     assert response.json() == new_molecule
#
#
# @pytest.mark.asyncio
# async def test_add_molecule_with_existing_id(async_client):
#     new_molecule = {
#         "id": 1,
#         "name": "Duplicate Water",
#         "smiles": "O"
#     }
#     response = await async_client.post("/molecules/add", json=new_molecule)
#     assert response.status_code == 400
#     assert response.json() == {"detail": "Molecule with this ID already exists."}
#
#
# @pytest.mark.asyncio
# async def test_add_molecule_with_existing_name(async_client):
#     new_molecule = {
#         "id": 11,
#         "name": "Water",
#         "smiles": "O"
#     }
#     response = await async_client.post("/molecules/add", json=new_molecule)
#     assert response.status_code == 400
#     assert response.json() == {"detail": "Molecule with this name already exists."}
#
#
# @pytest.mark.asyncio
# async def test_list_molecules(async_client):
#     response = await async_client.get("/molecules/")
#     assert response.status_code == 200
#     assert isinstance(response.json(), list)
#
#
# @pytest.mark.asyncio
# async def test_get_molecule_by_id(async_client):
#     response = await async_client.get("/molecules/1")
#     assert response.status_code == 200
#     assert response.json()["id"] == 1
#
#
# @pytest.mark.asyncio
# async def test_update_molecule_by_id(async_client):
#     updated_molecule = {
#         "id": 1,
#         "name": "Updated Water",
#         "smiles": "O"
#     }
#     response = await async_client.put("/molecules/1", json=updated_molecule)
#     assert response.status_code == 200
#     assert response.json() == updated_molecule
#
#
# @pytest.mark.asyncio
# async def test_update_nonexistent_molecule_by_id(async_client):
#     updated_molecule = {
#         "id": 999,
#         "name": "Nonexistent",
#         "smiles": "N"
#     }
#     response = await async_client.put("/molecules/999", json=updated_molecule)
#     assert response.status_code == 404
#     assert response.json() == {"detail": "Molecule not found."}
#
#
# @pytest.mark.asyncio
# async def test_delete_molecule_by_id(async_client):
#     response = await async_client.delete("/molecules/1")
#     assert response.status_code == 200
#     assert response.json()["id"] == 1
#
#
# @pytest.mark.asyncio
# async def test_delete_nonexistent_molecule_by_id(async_client):
#     response = await async_client.delete("/molecules/999")
#     assert response.status_code == 404
#     assert response.json() == {"detail": "Molecule not found."}
#
#
# @pytest.mark.asyncio
# async def test_search_molecules_by_substructure(async_client):
#     response = await async_client.get("/molecules/search/?mol=c1ccccc1")
#     assert response.status_code == 200
#     assert "matches" in response.json()
#     assert response.json()["matches"] == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
#
#
# @pytest.mark.asyncio
# async def test_search_with_empty_substructure(async_client):
#     response = await async_client.get("/molecules/search/?mol=")
#     assert response.status_code == 400
#     assert response.json() == {"detail": "Invalid substructure SMILES string. Must be a non-empty string."}
#
#
# @pytest.mark.asyncio
# async def test_upload_file_valid(async_client):
#     valid_smiles = "CCO\nCCN\nc1ccccc1\n"
#     response = await async_client.post(
#         "/molecules/upload/",
#         files={"file": ("valid_smiles.txt", valid_smiles)},
#     )
#     assert response.status_code == 201
#     assert "molecules" in response.json()
#     assert len(response.json()["molecules"]) == 3
#
#
# @pytest.mark.asyncio
# async def test_upload_file_invalid_type(async_client):
#     invalid_file = "Not a SMILES string"
#     response = await async_client.post(
#         "/molecules/upload/",
#         files={"file": ("invalid_file.csv", invalid_file)},
#     )
#     assert response.status_code == 400
#     assert response.json() == {"detail": "Invalid file type. Please upload a .txt file."}
#
#
# @pytest.mark.asyncio
# async def test_upload_file_no_valid_smiles(async_client):
#     invalid_smiles = "invalid_smiles\nanother_invalid_smile\n"
#     response = await async_client.post(
#         "/molecules/upload/",
#         files={"file": ("invalid_smiles.txt", invalid_smiles)},
#     )
#     assert response.status_code == 422
#     assert response.json() == {"detail": "File did not contain a single valid SMILES string."}
