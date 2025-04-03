#!/bin/bash

# Baue docker image
docker build -t bloodagent .

#teste docker image
docker run --rm -it bloodagent

# Baue singularity image
singularity build bloodagent.sif docker-daemon://bloodagent:latest

# test
singularity exec bloodagent.sif /app/bloodAGENT

# AUngetaggte Images löschen
docker image prune

