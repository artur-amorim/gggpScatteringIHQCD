CXX=g++
CXXFLAGS=-Iinclude -std=c++11 -pthread -O3 -c
LIBS=-larmadillo -lgsl -lblas -lpthread
LIBS_GFORTRAN = -L/usr/local/lib/gcc/9 -lgfortran

FC = gfortran

all: mkfolder tests

clean:
	rm -r bin

# test code

tests: bin/test_U1NNMode.exe bin/test_schrodinger.exe bin/test_SoftPomeron.exe bin/test_HardPomeron.exe bin/test_alphaQED.exe bin/test_PhotonScattering.exe bin/test_SigmaGammaGamma.exe

bin/test_SigmaGammaGamma.exe: bin/tests/test_SigmaGammaGamma.o bin/PhotonScattering.o bin/SigmaGammaGamma.o bin/Reggeon.o bin/Spectra.o bin/IHQCD.o bin/quadpack/dqags.o bin/quadpack/dqagse.o bin/quadpack/dqelg.o bin/quadpack/dqk21.o bin/quadpack/dqpsrt.o bin/quadpack/xerror.o bin/quadpack/xerrwv.o bin/quadpack/fdump.o bin/quadpack/j4save.o bin/quadpack/xerabt.o bin/quadpack/xerctl.o bin/quadpack/xerprt.o bin/quadpack/xersav.o bin/quadpack/xgetua.o bin/quadpack/s88fmt.o bin/quadpack/d1mach.o bin/quadpack/i1mach.o
	$(CXX) -o $@ $^ $(LIBS) $(LIBS_GFORTRAN)

bin/test_PhotonScattering.exe: bin/tests/test_PhotonScattering.o bin/PhotonScattering.o bin/Reggeon.o bin/Spectra.o
	$(CXX) -o $@ $^ $(LIBS)

bin/test_alphaQED.exe: bin/tests/test_alphaQED.o bin/alphaQED.o
	$(CXX) -o $@ $^ $(LIBS)

bin/test_HardPomeron.exe: bin/tests/test_HardPomeron.o bin/IHQCD.o bin/Reggeon.o bin/Kernel.o bin/HardPomeron.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/chebspec.o bin/schrodinger/numerov.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)

bin/test_SoftPomeron.exe: bin/tests/test_SoftPomeron.o bin/IHQCD.o bin/Reggeon.o bin/Kernel.o bin/SoftPomeron.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/chebspec.o bin/schrodinger/numerov.o bin/schrodinger/schrodinger.o
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