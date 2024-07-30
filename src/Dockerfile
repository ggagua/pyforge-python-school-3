FROM continuumio/miniconda3

WORKDIR /app

#Had to install pip version because conda did not work, it had problems with 3.12 Python version
RUN pip install rdkit

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
