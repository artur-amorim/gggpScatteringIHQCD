CXX=g++
CXXFLAGS=-Iinclude -std=c++11 -pthread -O3 -c
LIBS=-larmadillo -lgsl -lblas -lpthread
LIBS_GFORTRAN = -L/usr/local/lib/gcc/9 -lgfortran

FC = gfortran

all: mkfolder tests

clean:
	rm -r bin

# test code

tests: bin/test_U1NNMode.exe bin/test_schrodinger.exe bin/test_Kernel.exe bin/test_alphaQED.exe bin/test_PhotonScattering.exe

bin/test_PhotonScattering.exe: bin/tests/test_PhotonScattering.o bin/PhotonScattering.o bin/Reggeon.o bin/Spectra.o
	$(CXX) -o $@ $^ $(LIBS)

bin/test_alphaQED.exe: bin/tests/test_alphaQED.o bin/alphaQED.o
	$(CXX) -o $@ $^ $(LIBS)

bin/test_Kernel.exe: bin/tests/test_Kernel.o bin/IHQCD.o bin/Reggeon.o bin/Kernel.o bin/GluonKernel.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/chebspec.o bin/schrodinger/numerov.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)

bin/test_U1NNMode.exe: bin/tests/test_U1NNMode.o bin/U1NNMode.o bin/IHQCD.o bin/Fortran/colnew.o bin/Fortran/dgefa.o bin/Fortran/dgesl.o
	$(CXX) -o $@ $^ $(LIBS) $(LIBS_GFORTRAN)

bin/test_schrodinger.exe: bin/tests/test_schrodinger.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/chebspec.o bin/schrodinger/numerov.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)
# test files

bin/tests/%.o: tests/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

# src files

# Fortran
bin/Fortran/%.o: src/Fortran/%.f
	$(FC) $^ -std=legacy -c
	mv *.o bin/Fortran/

bin/quadpack/%.o: src/quadpack/%.f
	$(FC) $^ -std=legacy -c
	mv *.o bin/quadpack/

# C++ files
bin/%.o: src/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

mkfolder: bin bin/schrodinger bin/tests bin/Fortran bin/quadpack bin/physics

bin/Fortran:
	mkdir $@

bin/tests:
	mkdir $@

bin/quadpack:
	mkdir $@

bin/physics:
	mkdir $@

bin/schrodinger:
	mkdir $@
bin:
	mkdir $@