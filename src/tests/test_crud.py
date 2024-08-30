# import pytest
# from src import crud, schemas
# from sqlalchemy.ext.asyncio import AsyncSession
#
# @pytest.mark.asyncio
# async def test_create_molecule(db: AsyncSession):
#     molecule_data = schemas.MoleculeCreate(name="Test Molecule", smiles="CCO")
#     result = await crud.create_molecule(db, molecule_data)
#     assert result.name == molecule_data.name
#     assert result.smiles == molecule_data.smiles