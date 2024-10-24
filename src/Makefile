# Simple Makefile used to compile a collection of C++ API functions
# contained in the specified $(LIBNAME).cpp as a "shared object" library.

LIBNAME=ParticleDynamicsAPI
TESTNAME=LineLoadTest
EXONAME=example_space_frame.exo

CC=g++
FLAGS=-O2 -Wall -pedantic -shared -fPIC -lstdc++ -std=c++11

.PHONY: clean debug test

all: $(LIBNAME).so

test: $(LIBNAME).so $(TESTNAME).py $(EXONAME) ParticleDynamics.py 
	python3 $(TESTNAME).py --file $(EXONAME)

debug: $(LIBNAME).cpp
	$(CC) $(FLAGS) -DDEBUG_ENABLE -o $(LIBNAME).so $(LIBNAME).cpp

$(LIBNAME).so: $(LIBNAME).cpp
	$(CC) $(FLAGS) -o $(LIBNAME).so $(LIBNAME).cpp

clean:
	rm -rf $(LIBNAME).so __pycache__ $(TESTNAME)_PVD $(TESTNAME)_PVD.pvd particles.exo
