FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md LICENSE ./
COPY frap_toolbox_py ./frap_toolbox_py

RUN python -m pip install --upgrade pip \
    && python -m pip install ".[app-full,cloud]"

ENTRYPOINT ["frap-toolbox-run"]
