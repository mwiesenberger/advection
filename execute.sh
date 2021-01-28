#!/bin/bash

# Script to run a simulation

: ${FELTOR_PATH:="../feltor"}
# If feltor is not here then change the FELTOR_PATH enviromnent variable
# export FELTOR_PATH="path/to/feltor"

echo "Compiling the source code ... "
make -C $FELTOR_PATH/src/lamb_dipole shu_b
echo "... Done"

echo "$FELTOR_PATH/src/lamb_dipole/shu_b $1 $2"
$FELTOR_PATH/src/lamb_dipole/shu_b $1 $2
