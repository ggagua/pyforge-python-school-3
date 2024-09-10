from src.database import Base
from sqlalchemy.orm import Mapped, mapped_column


class Molecule(Base):
    __tablename__ = 'molecules'

    id: Mapped[int] = mapped_column(primary_key=True, index=True)
    name: Mapped[str] = mapped_column(unique=True, index=True)
    smiles: Mapped[str] = mapped_column(unique=True, index=True)
