#!/bin/bash
set -xuo pipefail

for TAG in $(git tag -l)
  do
    asv run --skip-existing-successful "$@" "${TAG}^..${TAG}"
  done
