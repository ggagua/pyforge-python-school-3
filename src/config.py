from pydantic import PostgresDsn
from pydantic_settings import BaseSettings
from dotenv import load_dotenv

load_dotenv()


# Inspired from https://github.com/fastapi/full-stack-fastapi-template/blob/master/backend/app/core/config.py
class Settings(BaseSettings):
    PROJECT_NAME: str
    SENTRY_DSN: str | None = None
    DB_SERVER: str = "localhost"  # Use your default value or adjust as needed in your .env file
    DB_PORT: int = 5433
    DB_USER: str
    DB_PASSWORD: str
    DB_NAME: str

    @property
    def SQLALCHEMY_DATABASE_URI(self) -> str:
        return f"postgresql+asyncpg://{self.DB_USER}:{self.DB_PASSWORD}@{self.DB_SERVER}:{self.DB_PORT}/{self.DB_NAME}"

    class Config:
        env_file = ".env"


settings = Settings()
