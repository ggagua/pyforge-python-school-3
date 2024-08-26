# FastAPI Imports
from src.database import Base
from pydantic import BaseModel
from sqlalchemy.orm import Mapped, mapped_column



class Molecule(Base):
    __tablename__ = 'molecules'

    id: Mapped[int] = mapped_column(primary_key=True, index=True)
    name: Mapped[str] = mapped_column(unique=True, index=True)
    smiles: Mapped[str] = mapped_column(unique=True, index=True)


class MoleculeCreate(BaseModel):
    name: str
    smiles: str


class MoleculeInDB(BaseModel):
    id: int
    name: str
    smiles: str

    class Config:
        orm_mode = True


class MoleculeUpdate(BaseModel):
    name: str
    smiles: str


class MoleculeResponse(BaseModel):
    id: int
    name: str
    smiles: str

    class Config:
        orm_mode = True
