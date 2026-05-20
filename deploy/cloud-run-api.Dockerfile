FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY frap_toolbox_py ./frap_toolbox_py

RUN python -m pip install --upgrade pip \
    && python -m pip install ".[cloud-api,cloud-gcp]"

CMD ["uvicorn", "frap_toolbox_py.cloud.backend:app", "--host", "0.0.0.0", "--port", "8080"]
