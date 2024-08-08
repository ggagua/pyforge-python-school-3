import pytest
from fastapi.testclient import TestClient
from .main import substructure_search, app

client = TestClient(app)


# Parameterized Testing for substructure_search with valid inputs
@pytest.mark.parametrize(
    "mols, substructure, expected",
    [
        (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1", ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),
        (["CCO", "CCC", "CCN"], "C", ["CCO", "CCC", "CCN"]),
        (["CCO", "CCC", "CCN"], "N", ["CCN"]),
        (["CCO", "CCC", "CCN"], "O", ["CCO"]),
        (["CCO", "CCC"], "CC", ["CCO", "CCC"]),
    ]
)
def test_substructure_search_param(mols, substructure, expected):
    assert substructure_search(mols, substructure) == expected


# Testing for invalid substructure SMILES string
def test_substructure_search_invalid_substructure():
    with pytest.raises(ValueError, match="Unknown substructure SMILES string"):
        substructure_search(["CCO", "c1ccccc1"], "invalid_smiles")


# Testing for invalid molecule SMILES string
def test_substructure_search_invalid_molecule_smiles():
    mols = ["CCO", "CCC", "invalid_smiles"]
    substructure = "CC"
    with pytest.raises(ValueError, match="Invalid molecule SMILES string: invalid_smiles"):
        substructure_search(mols, substructure)


# Testing for error processing molecule
def test_substructure_search_error_processing_molecule():
    mols = ["CCO", "CCC", "invalid_smiles"]
    substructure = "CC"
    with pytest.raises(ValueError,
                       match="Error processing molecule invalid_smiles: Invalid molecule SMILES string: invalid_smiles"):
        substructure_search(mols, substructure)


# Testing FastAPI endpoints
def test_get_server():
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()


def test_add_molecule():
    new_molecule = {
        "id": 10,
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }
    response = client.post("/molecules/add", json=new_molecule)
    assert response.status_code == 201
    assert response.json() == new_molecule


def test_add_molecule_with_existing_id():
    new_molecule = {
        "id": 1,
        "name": "Duplicate Water",
        "smiles": "O"
    }
    response = client.post("/molecules/add", json=new_molecule)
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule with this ID already exists."}


def test_add_molecule_with_existing_name():
    new_molecule = {
        "id": 11,
        "name": "Water",
        "smiles": "O"
    }
    response = client.post("/molecules/add", json=new_molecule)
    assert response.status_code == 400
    assert response.json() == {"detail": "Molecule with this name already exists."}


def test_list_molecules():
    response = client.get("/molecules/")
    assert response.status_code == 200
    assert isinstance(response.json(), list)


def test_get_molecule_by_id():
    response = client.get("/molecules/1")
    assert response.status_code == 200
    assert response.json()["id"] == 1


def test_update_molecule_by_id():
    updated_molecule = {
        "id": 1,
        "name": "Updated Water",
        "smiles": "O"
    }
    response = client.put("/molecules/1", json=updated_molecule)
    assert response.status_code == 200
    assert response.json() == updated_molecule


def test_update_nonexistent_molecule_by_id():
    updated_molecule = {
        "id": 999,
        "name": "Nonexistent",
        "smiles": "N"
    }
    response = client.put("/molecules/999", json=updated_molecule)
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found."}


def test_delete_molecule_by_id():
    response = client.delete("/molecules/1")
    assert response.status_code == 200
    assert response.json()["id"] == 1


def test_delete_nonexistent_molecule_by_id():
    response = client.delete("/molecules/999")
    assert response.status_code == 404
    assert response.json() == {"detail": "Molecule not found."}


def test_search_molecules_by_substructure():
    response = client.get("/molecules/search/?mol=c1ccccc1")
    assert response.status_code == 200
    assert "matches" in response.json()
    assert response.json()["matches"] == ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]


def test_search_with_empty_substructure():
    response = client.get("/molecules/search/?mol=")
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid substructure SMILES string. Must be a non-empty string."}


# Test for uploading a valid text file containing SMILES strings
def test_upload_file_valid():
    # Creating a valid text file with SMILES strings
    valid_smiles = "CCO\nCCN\nc1ccccc1\n"
    response = client.post(
        "/molecules/upload/",
        files={"file": ("valid_smiles.txt", valid_smiles)},
    )
    assert response.status_code == 201
    assert "molecules" in response.json()
    assert len(response.json()["molecules"]) == 3  # Assuming all are valid


# Test for uploading an invalid file type
def test_upload_file_invalid_type():
    # Create an invalid file type (e.g., a .csv file)
    invalid_file = "Not a SMILES string"
    response = client.post(
        "/molecules/upload/",
        files={"file": ("invalid_file.csv", invalid_file)},
    )
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid file type. Please upload a .txt file."}


# Test for uploading a text file with no valid SMILES strings
def test_upload_file_no_valid_smiles():
    # Creating a text file with invalid strings
    invalid_smiles = "invalid_smiles\nanother_invalid_smile\n"
    response = client.post(
        "/molecules/upload/",
        files={"file": ("invalid_smiles.txt", invalid_smiles)},
    )
    assert response.status_code == 422
    assert response.json() == {"detail": "File did not contain a single valid SMILES string."}


# Test for getting an image by its name
def test_get_image_valid():
    # Assume the upload was successful and the image was generated
    valid_smiles = "CCO\nCCN\nc1ccccc1\n"
    upload_response = client.post(
        "/molecules/upload/",
        files={"file": ("valid_smiles.txt", valid_smiles)},
    )

    # Check that we got images created
    assert upload_response.status_code == 201
    image_name = upload_response.json()["molecules"][0]["image_name"]

    # Getting the image
    response = client.get(f"/molecules/images/{image_name}")
    assert response.status_code == 200
    assert response.headers["content-type"] == "image/png"


# Test for getting an image that does not exist
def test_get_image_not_found():
    response = client.get("/molecules/images/non_existent_image.png")
    assert response.status_code == 404
    assert response.json() == {"detail": "Image not found"}
