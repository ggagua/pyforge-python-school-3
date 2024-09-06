from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select
from sqlalchemy.exc import IntegrityError, SQLAlchemyError
from fastapi import HTTPException, status
from rdkit import Chem
from typing import List, Optional, AsyncIterator
from src import models, schemas


async def get_molecule(db: AsyncSession, molecule_id: int) -> Optional[models.Molecule]:
    """
    Retrieve a single molecule by its ID.
    """
    try:
        result = await db.execute(select(models.Molecule).filter(models.Molecule.id == molecule_id))
        return result.scalar_one_or_none()
    except SQLAlchemyError as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail=f"Error fetching molecule: {str(e)}")


async def get_molecules(db: AsyncSession, skip: int = 0, limit: int = 100) -> AsyncIterator[models.Molecule]:
    """
    Retrieve a list of molecules with optional pagination as an asynchronous generator.
    """
    try:
        result = await db.execute(select(models.Molecule).offset(skip).limit(limit))
        molecules = result.scalars().all()

        for molecule in molecules:
            yield molecule
    except SQLAlchemyError as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail=f"Error fetching molecules: {str(e)}")


async def create_molecule(db: AsyncSession, molecule: schemas.MoleculeCreate) -> models.Molecule:
    """
    Create a new molecule in the database.
    """
    db_molecule = models.Molecule(name=molecule.name, smiles=molecule.smiles)
    db.add(db_molecule)

    try:
        await db.commit()
        await db.refresh(db_molecule)
    except IntegrityError:
        await db.rollback()
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Molecule with this name or SMILES already exists."
        )
    except SQLAlchemyError as e:
        await db.rollback()
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail=f"Error creating molecule: {str(e)}")

    return db_molecule


async def update_molecule(db: AsyncSession, molecule_id: int, molecule: schemas.MoleculeUpdate) -> Optional[
    models.Molecule]:
    """
    Update an existing molecule by its ID.
    """
    try:
        result = await db.execute(select(models.Molecule).filter(models.Molecule.id == molecule_id))
        db_molecule = result.scalar_one_or_none()

        if db_molecule:
            db_molecule.name = molecule.name
            db_molecule.smiles = molecule.smiles

            try:
                await db.commit()
                await db.refresh(db_molecule)
            except IntegrityError:
                await db.rollback()
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail="Updated molecule would violate unique constraints."
                )
            except SQLAlchemyError as e:
                await db.rollback()
                raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                                    detail=f"Error updating molecule: {str(e)}")

            return db_molecule

        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")
    except SQLAlchemyError as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail=f"Error fetching molecule: {str(e)}")


async def delete_molecule(db: AsyncSession, molecule_id: int) -> Optional[models.Molecule]:
    """
    Delete a molecule by its ID.
    """
    try:
        result = await db.execute(select(models.Molecule).filter(models.Molecule.id == molecule_id))
        db_molecule = result.scalar_one_or_none()

        if db_molecule:
            try:
                await db.delete(db_molecule)
                await db.commit()
                return db_molecule
            except SQLAlchemyError as e:
                await db.rollback()
                raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                                    detail=f"Error deleting molecule: {str(e)}")

        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")
    except SQLAlchemyError as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail=f"Error fetching molecule: {str(e)}")


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
