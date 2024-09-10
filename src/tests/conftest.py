# import pytest
# from httpx import AsyncClient
# from sqlalchemy.ext.asyncio import AsyncSession, create_async_engine
# from sqlalchemy.ext.asyncio import async_sessionmaker
# from sqlalchemy.orm import sessionmaker
# from src.main import app
# from src.config import settings
# from src.database import init_db, drop_test_db
#
# # Async engine and session setup
# engine = create_async_engine(settings.SQLALCHEMY_DATABASE_URI, echo=True)
# TestingSessionLocal = async_sessionmaker(bind=engine, class_=AsyncSession)
#
# @pytest.fixture(scope="session", autouse=True)
# async def setup_and_teardown():
#     # Create test database and tables
#     async with engine.begin() as conn:
#         await drop_test_db(engine)
#         await init_db()
#     yield
#     # Drop test database and tables after tests
#     async with engine.begin() as conn:
#         await drop_test_db(engine)
#
# @pytest.fixture
# async def client():
#     async with AsyncClient(
#             transport=ASGITransport(app=app), base_url="http://test"
#     ) as ac:
#         yield ac