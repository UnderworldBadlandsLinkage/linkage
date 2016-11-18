# UnderworldBadlandsLinkage

Run combined Underworld and Badlands models

## Getting started

Grab the <todo> container from Docker Hub. Start it up, and you should get a Jupyter Notebook with a selection of demonstration models.

## Docker containers

We provide three Docker containers:

docker-base: The third-part dependencies required by Badlands and Underworld. This does not change very often.
docker-models: Builds on docker-base and installs the Badlands and Underworld models. This changes fairly often at the moment.
docker-demo: A container that add some demonstration Jupyter Notebooks and preconfigures an MPI cluster for you. **Unless you're developing Badlands or Underworld, you should be using this one.**

We split the containers to reduce download size when we update the models.

