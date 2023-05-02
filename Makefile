device=cpu
FELTOR_PATH=../feltor

#configure machine
include $(FELTOR_PATH)/config/default.mk
include $(FELTOR_PATH)/config/version.mk
include $(FELTOR_PATH)/config/*.mk
include $(FELTOR_PATH)/config/devices/devices.mk

INCLUDE+=-I$(FELTOR_PATH)/inc/

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
