import json
from unittest.mock import AsyncMock, create_autospec
from fastapi import HTTPException, status
import pytest
from fastapi.testclient import TestClient
from src import crud
from src.crud import substructure_search
from src.schemas import MoleculeCreate, MoleculeUpdate, MoleculeSearch
from src.main import app
from unittest.mock import MagicMock, patch
from rdkit import Chem
from unittest import mock
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.engine.result import ScalarResult

client = TestClient(app)


# Decided to go with mocking database instead of connecting to it, since I've had infinite amount of problems with creating
# test Postgres Database, migrating to it and testing it. Did not want to use sqlite as well...
# Will improve these tests later as well.
@pytest.fixture
def mock_create_molecule(monkeypatch):
    async def mock_post(db, molecule: MoleculeCreate):
        return {"id": 1, "name": molecule.name, "smiles": molecule.smiles}

    monkeypatch.setattr("src.crud.create_molecule", AsyncMock(side_effect=mock_post))


@pytest.mark.asyncio
async def test_add_molecule(mock_create_molecule):
    # Mock request data
    response = client.post("/molecules/add", json={"name": "something", "smiles": "CCCCCCHC"})

    assert response.status_code == 200
    assert response.json() == {"id": 1, "name": "something", "smiles": "CCCCCCHC"}


@pytest.fixture
def mock_create_molecule_conflict(monkeypatch):
    async def mock_post(db, molecule: MoleculeCreate):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Molecule with this name or SMILES already exists."
        )

    monkeypatch.setattr("src.crud.create_molecule", AsyncMock(side_effect=mock_post))


@pytest.mark.asyncio
async def test_create_molecule_conflict(mock_create_molecule_conflict):
    response = client.post("/molecules/add", json={"name": "something", "smiles": "CCCC"})
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule with this name or SMILES already exists."}


@pytest.fixture
def mock_get_molecule(monkeypatch):
    async def mock_get(db, molecule_id: int):
        return {"id": molecule_id, "name": "something", "smiles": "CCCCCCHC"}

    # Apply the monkeypatch for the `get_molecule` function
    monkeypatch.setattr("src.crud.get_molecule", AsyncMock(side_effect=mock_get))


@pytest.mark.asyncio
async def test_get_molecule(mock_get_molecule):
    molecule_id = 1
    response = client.get(f"/molecules/{molecule_id}")

    assert response.status_code == 200
    assert response.json() == {"id": 1, "name": "something", "smiles": "CCCCCCHC"}


@pytest.fixture
def mock_update_molecule(monkeypatch):
    async def mock_update(db, molecule_id: int, molecule: MoleculeUpdate):
        return {"id": molecule_id, "name": molecule.name, "smiles": molecule.smiles}

    monkeypatch.setattr("src.crud.update_molecule", AsyncMock(side_effect=mock_update))


@pytest.mark.asyncio
async def test_update_molecule(mock_update_molecule):
    molecule_id = 1
    updated_data = {"name": "updated_name", "smiles": "CCCC"}

    response = client.put(f"/molecules/{molecule_id}", json=updated_data)

    assert response.status_code == 200

    assert response.json() == {"id": molecule_id, "name": "updated_name", "smiles": "CCCC"}


@pytest.fixture
def mock_delete_molecule(monkeypatch):
    async def mock_delete(db, molecule_id: int):
        return {"id": molecule_id, "name": 'updated_molecule', "smiles": 'CCCCH'}

    monkeypatch.setattr("src.crud.delete_molecule", AsyncMock(side_effect=mock_delete))


@pytest.mark.asyncio
async def test_invalid_delete_molecule(mock_delete_molecule):
    molecule_id = 1
    response = client.delete(f"/molecules/{molecule_id}")
    assert response.status_code == 200
    assert response.json() == {"id": molecule_id, "name": "updated_molecule", "smiles": "CCCCH"}

# WILL FIGURE OUT SUBSTRUCTURE SEARCH LATER, HAD SOME PROBLEMS WITH WRITING IT.
