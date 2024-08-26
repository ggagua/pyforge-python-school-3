from fastapi import FastAPI
from src.routes.molecules_router import router as molecules_router
from src.routes.files_router import router as file_router
from os import getenv

app = FastAPI()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


app.include_router(molecules_router, prefix="/molecules", tags=["Molecules"])
app.include_router(file_router, prefix="/files", tags=["Files"])
