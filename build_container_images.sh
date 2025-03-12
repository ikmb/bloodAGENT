#!/bin/bash

# Baue docker image
docker build -t bloodagent .

#teste docker image
docker run --rm -it bloodagent

# Baue singularity image
singularity build bloodagent.sif docker-daemon://bloodagent:latest



# Alles unter Docker löschen (es wird gefragt)
docker system prune -a --volumes

