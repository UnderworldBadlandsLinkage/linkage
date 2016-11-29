#!/bin/bash

# Modify this script to suit your operating environment

# Path to the base Underworld directory
UNDERWORLD_PATH=/Users/something/underworld2

# Path to the 'libUtils' directory under pyBadlands
LIBUTILS_PATH=/Users/something/pyBadlands/pyBadlands/libUtils

# before you start
# Xvfb :0 -screen 0 1600x1200x16&

# if you're using visualisation
#export DISPLAY=:0

export PYTHONPATH=$PYTHONPATH:$UNDERWORLD_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBUTILS_PATH

# The simplest way to run is just 'python <script name>'
cd dem
python dem.py

# If you're using a virtualenv, activate it before running the script
#bash -c '. ../../env/bin/activate && cd dem && python dem.py'

# You can run in MPI with mpirun
# mpirun -n 4 'cd dem && python dem.py'

# And of course, you can use a virtualenv with MPI:
# mpirun -n 4 bash -c '. ../../env/bin/activate && cd dem && python dem.py'
