FROM python:3.5
RUN apt-get update && apt-get install -yq --no-install-recommends \
    vim \
    postgresql \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN mkdir /mc
WORKDIR /mc/dj
ADD requirements.txt /requirements.txt
RUN pip install -r /requirements.txt
