#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RDKIT_CONFORMER_GENERATOR_DIR=$(dirname $DIR)

IMAGE_NAME="rdkit_conformer_generator_test:latest"
if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
  docker build . -t $IMAGE_NAME
fi

docker run \
  -v $RDKIT_CONFORMER_GENERATOR_DIR:/opt/rdkit_conformer_generator \
  --rm \
  -it \
  $IMAGE_NAME \
  /bin/bash -c "cd /opt; python -m rdkit_conformer_generator.tests.test_rdkit_conformer_generator $@"
