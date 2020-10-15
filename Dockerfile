FROM informaticsmatters/rdkit-python3-debian:latest
WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt
RUN pip install mordred
COPY /app .

CMD ["python", "main.py"]