import pytest
from sqlalchemy.ext.asyncio import AsyncSession
from src.database import get_db, engine, Base

#???????
#Attemptin but not working
@pytest.fixture(scope="function")
async def db() -> AsyncSession:
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)

    async with get_db() as session:
        yield session

    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.drop_all)