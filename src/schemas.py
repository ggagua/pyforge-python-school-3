from pydantic import BaseModel


# Will add some later
class MoleculeBase(BaseModel):
    name: str
    smiles: str


class MoleculeCreate(MoleculeBase):
    pass


class Molecule(MoleculeBase):
    id: int

    class Config:
        from_attributes = True
        json_schema_extra = {
            "examples": [
                {
                    "id": 10,
                    "name": "Caffeine",
                    "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
                },
            ]
        }


class MoleculeUpdate(MoleculeBase):
    pass


class MoleculeSearch(BaseModel):
    substructure: str
