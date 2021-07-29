#!/bin/bash

: ${FELTOR_PATH:="../feltor"}
# If feltor is not here then change the FELTOR_PATH enviromnent variable
# export FELTOR_PATH="path/to/feltor"

make -C $FELTOR_PATH/src/navier_stokes/ navier_stokes

$FELTOR_PATH/src/navier_stokes/navier_stokes $@
