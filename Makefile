: ${FELTOR_PATH:="../feltor"}

CC=g++ #C++ compiler
CFLAGS=-Wall -std=c++14 -mavx -mfma -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP
OPT=-O2 -fopenmp # optimization flag
# it is O2 and not O3 because g++-7 up to g++-8.0 have a bug with fma in -O3,
# fixed in g++-8.1

#external libraries
INCLUDE = -I$(HOME)/include    # cusp, thrust, jsoncpp and feltor/inc/dg links
LIBS=-lnetcdf -lhdf5 -lhdf5_hl  # netcdf library for file output
JSONLIB=-L$(HOME)/include/json/../../src/lib_json -ljsoncpp # json library for input parameters

all: continuity navier_stokes plasma

continuity: continuity.cpp common.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -g

navier_stokes: navier_stokes.cpp common.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -g

plasma: plasma.cpp common.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -g

.PHONY: clean

clean:
	rm -f continuity navier_stokes plasma
