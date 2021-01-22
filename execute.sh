#!/bin/bash

# Script to run a simulation

# Use
# export FELTOR_PATH="path/to/feltor"
# in the environment to set the feltor path to
# a different value than the default one

: ${FELTOR_PATH:="../feltor"}

echo "Compiling the source code ... "
make $FELTOR_PATH/src/lamb_dipole/Makefile
echo "... Done"

echo "$FELTOR_PATH/src/lamb_dipole/shu_b $1 $2"
$FELTOR_PATH/src/lamb_dipole/shu_b $1 $2
