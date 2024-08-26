from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select
from . import models, schemas
from rdkit import Chem
from typing import List, Optional

async def get_molecule(db: AsyncSession, molecule_id: int) -> Optional[models.Molecule]:
    """
    Retrieve a single molecule by its ID.

    :param db: The database session.
    :param molecule_id: The ID of the molecule to retrieve.
    :return: The molecule with the specified ID, or None if not found.
    """
    result = await db.execute(select(models.Molecule).filter(models.Molecule.id == molecule_id))
    return result.scalar_one_or_none()

async def get_molecules(db: AsyncSession, skip: int = 0) -> List[models.Molecule]:
    """
    Retrieve a list of molecules with optional pagination.

    :param db: The database session.
    :param skip: Number of records to skip (for pagination).
    :return: A list of molecules.
    """
    result = await db.execute(select(models.Molecule).offset(skip))
    return result.scalars().all()

async def create_molecule(db: AsyncSession, molecule: schemas.MoleculeCreate) -> models.Molecule:
    """
    Create a new molecule in the database.

    :param db: The database session.
    :param molecule: The data to create the new molecule.
    :return: The created molecule.
    """
    db_molecule = models.Molecule(name=molecule.name, smiles=molecule.smiles)
    db.add(db_molecule)
    await db.commit()
    await db.refresh(db_molecule)
    return db_molecule

async def update_molecule(db: AsyncSession, molecule_id: int, molecule: schemas.MoleculeUpdate) -> Optional[models.Molecule]:
    """
    Update an existing molecule by its ID.

    :param db: The database session.
    :param molecule_id: The ID of the molecule to update.
    :param molecule: The updated data for the molecule.
    :return: The updated molecule, or None if not found.
    """
    result = await db.execute(select(models.Molecule).filter(models.Molecule.id == molecule_id))
    db_molecule = result.scalar_one_or_none()
    if db_molecule:
        db_molecule.name = molecule.name
        db_molecule.smiles = molecule.smiles
        await db.commit()
        await db.refresh(db_molecule)
        return db_molecule
    return None

async def delete_molecule(db: AsyncSession, molecule_id: int) -> Optional[models.Molecule]:
    """
    Delete a molecule by its ID.

    :param db: The database session.
    :param molecule_id: The ID of the molecule to delete.
    :return: The deleted molecule, or None if not found.
    """
    result = await db.execute(select(models.Molecule).filter(models.Molecule.id == molecule_id))
    db_molecule = result.scalar_one_or_none()
    if db_molecule:
        await db.delete(db_molecule)
        await db.commit()
        return db_molecule
    return None


async def substructure_search(db: AsyncSession, mol: str) -> List[str]:
    """
    Search for molecules that contain a given substructure.

    :param db: The database session.
    :param mol: The SMILES string representing the substructure.
    :return: A list of SMILES strings for molecules containing the substructure.
    """
    # Convert substructure SMILES to RDKit Mol object
    substructure = Chem.MolFromSmiles(mol)
    if substructure is None:
        raise ValueError("Invalid substructure SMILES string")

    result = await db.execute(select(models.Molecule))
    molecules = result.scalars().all()

    # Prepare to collect matches
    matches = []
    for molecule in molecules:
        try:
            # Convert molecule SMILES to RDKit Mol object
            molecule_smiles = molecule.smiles
            molecule_rdkit = Chem.MolFromSmiles(molecule_smiles)

            if molecule_rdkit is None:
                continue  # Skip invalid molecule SMILES strings

            # Perform substructure match
            if molecule_rdkit.HasSubstructMatch(substructure):
                matches.append(molecule_smiles)
        except Exception as e:
            raise ValueError(f"Error processing molecule with SMILES {molecule_smiles}: {e}")

    return matches