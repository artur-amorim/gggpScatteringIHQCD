# Specify the minimum CMAKE version required
cmake_minimum_required(VERSION 2.8.0)

set(CMAKE_Fortran_COMPILER gfortran)
ENABLE_LANGUAGE(Fortran)
set(CMAKE_Fortran_FLAGS "-std=legacy")

if(UNIX AND NOT APPLE)
        set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -pthread")
        FIND_PACKAGE(Boost 1.53 REQUIRED)
        include_directories(${Boost_INCLUDE_DIRS} "include")
endif()

if(APPLE)
        set(CMAKE_CXX_FLAGS "-std=c++11 -O3")   
        set(OS_INCLUDE_DIRS "/usr/include" "/usr/local/include")
        set(OS_LIB_DIRS "/usr/lib" "/usr/local/lib")
        include_directories(${OS_INCLUDE_DIRS} "include")
        link_directories(${OS_LIB_DIRS})
endif()

make_directory(execs)

set(EXECUTABLE_OUTPUT_PATH "execs")

# Variable SCHRODINGER to be used when schrodinger src code is needed
set(SCHRODINGER src/schrodinger/common.cpp src/schrodinger/solvspec.cpp src/schrodinger/chebspec.cpp src/schrodinger/numerov.cpp src/schrodinger/schrodinger.cpp)
#add_library(schrodinger_lib ${SCHRODINGER})

# QUADPACK
set(QUADPACK src/quadpack/d1mach.f src/quadpack/dqags.f src/quadpack/dqagse.f src/quadpack/dqelg.f src/quadpack/dqk21.f src/quadpack/dqpsrt.f src/quadpack/fdump.f src/quadpack/i1mach.f src/quadpack/j4save.f src/quadpack/s88fmt.f src/quadpack/xerabt.f src/quadpack/xerctl.f src/quadpack/xerprt.f src/quadpack/xerror.f src/quadpack/xerrwv.f src/quadpack/xersav.f src/quadpack/xgetua.f)
#add_library(quadpack_lib ${QUADPACK})

# COLNEW
set(COLNEW src/Fortran/colnew.f src/Fortran/dgefa.f src/Fortran/dgesl.f)
#add_library(colnew_lib ${COLNEW})

# Variable YangMills
set(IHQCD src/IHQCD.cpp)

# Variable HVQCD
set(HVQCD src/background.cpp src/HolographicVQCD.cpp)

# U1NNMode
set(U1NNMODE src/U1NNMode.cpp)

# Spectra and Reggeon classes
set(SPECTRA src/Spectra.cpp src/Reggeon.cpp src/Kernel.cpp src/HardPomeron.cpp src/SoftPomeron.cpp)

# ProcessObservable and its derived classes
set(PROCESS src/Process.cpp src/DeepInelasticScattering.cpp src/Sigma.cpp src/F2.cpp src/FL.cpp src/F2Photon.cpp src/SigmaGammaGamma.cpp src/SigmaGammaProton.cpp src/SigmaProtonProton.cpp src/alphaQED.cpp)

# HQCDP class
set(HQCDP src/HQCDP.cpp)

add_library(ihqcd_lib ${IHQCD} ${SPECTRA} ${PROCESS} ${U1NNMODE} ${HQCDP} ${SCHRODINGER} ${COLNEW} ${QUADPACK})

add_executable(test_alphaQED.exe tests/test_alphaQED.cpp)
target_link_libraries(test_alphaQED.exe ihqcd_lib)

add_executable(test_F2Photon.exe tests/test_F2Photon.cpp)
target_link_libraries(test_F2Photon.exe ihqcd_lib blas armadillo)

add_executable(test_HardPomeron.exe tests/test_HardPomeron.cpp)
target_link_libraries(test_HardPomeron.exe ihqcd_lib blas armadillo)

add_executable(test_SoftPomeron.exe tests/test_SoftPomeron.cpp)
target_link_libraries(test_SoftPomeron.exe ihqcd_lib blas armadillo)

add_executable(test_HQCDP_HardPomeron.exe tests/test_HQCDP_HardPomeron.cpp)
target_link_libraries(test_HQCDP_HardPomeron.exe ihqcd_lib blas armadillo)

add_executable(test_HQCDP_SoftPomeron.exe tests/test_HQCDP_SoftPomeron.cpp)
target_link_libraries(test_HQCDP_SoftPomeron.exe ihqcd_lib blas armadillo)

add_executable(test_schrodinger.exe tests/test_schrodinger.cpp)
target_link_libraries(test_schrodinger.exe ihqcd_lib blas armadillo)

add_executable(test_SigmaGammaGamma.exe tests/test_SigmaGammaGamma.cpp)
target_link_libraries(test_SigmaGammaGamma.exe ihqcd_lib blas armadillo)

add_executable(test_SigmaGammaProton.exe tests/test_SigmaGammaProton.cpp)
target_link_libraries(test_SigmaGammaProton.exe ihqcd_lib blas armadillo)

add_executable(test_U1NNMode.exe tests/test_U1NNMode.cpp)
target_link_libraries(test_U1NNMode.exe ihqcd_lib blas)

add_executable(fitF2Photon_HardPomeron_4Reggeons.exe physics/fitF2Photon_HardPomeron_4Reggeons.cpp)
target_link_libraries(fitF2Photon_HardPomeron_4Reggeons.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaGammaGamma_HardPomeron_2Reggeons.exe physics/fitSigmaGammaGamma_HardPomeron_2Reggeons.cpp)
target_link_libraries(fitSigmaGammaGamma_HardPomeron_2Reggeons.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaGammaGamma_HardPomeron_4Reggeons.exe physics/fitSigmaGammaGamma_HardPomeron_4Reggeons.cpp)
target_link_libraries(fitSigmaGammaGamma_HardPomeron_4Reggeons.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaGammaGamma_SoftPomeron.exe physics/fitSigmaGammaGamma_SoftPomeron.cpp)
target_link_libraries(fitSigmaGammaGamma_SoftPomeron.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaGammaGammaF2Photon_HardPomeron_4Reggeons.exe physics/fitSigmaGammaGammaF2Photon_HardPomeron_4Reggeons.cpp)
target_link_libraries(fitSigmaGammaGammaF2Photon_HardPomeron_4Reggeons.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaGammaProton_HardPomeron.exe physics/fitSigmaGammaProton_HardPomeron.cpp)
target_link_libraries(fitSigmaGammaProton_HardPomeron.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaGammaProton_SoftPomeron.exe physics/fitSigmaGammaProton_SoftPomeron.cpp)
target_link_libraries(fitSigmaGammaProton_SoftPomeron.exe ihqcd_lib blas armadillo)


add_executable(fitDIS_SigmaGammaProton.exe physics/fitDIS_SigmaGammaProton.cpp)
target_link_libraries(fitDIS_SigmaGammaProton.exe ihqcd_lib blas armadillo)

add_executable(predictDISStructureFunctions.exe physics/predictDISStructureFunctions.cpp)
target_link_libraries(predictDISStructureFunctions.exe ihqcd_lib blas armadillo)

add_executable(predictSigmaGammaGamma_HardPomeron_4Reggeons.exe physics/predictSigmaGammaGamma_HardPomeron_4Reggeons.cpp)
target_link_libraries(predictSigmaGammaGamma_HardPomeron_4Reggeons.exe ihqcd_lib blas armadillo)

add_executable(predictSigmaGammaGamma_SoftPomeron.exe physics/predictSigmaGammaGamma_SoftPomeron.cpp)
target_link_libraries(predictSigmaGammaGamma_SoftPomeron.exe ihqcd_lib blas armadillo)

add_executable(predictF2Photon.exe physics/predictF2Photon.cpp)
target_link_libraries(predictF2Photon.exe ihqcd_lib blas armadillo)

add_executable(predictSigmaGammaProton_HardPomeron.exe physics/predictSigmaGammaProton_HardPomeron.cpp)
target_link_libraries(predictSigmaGammaProton_HardPomeron.exe ihqcd_lib blas armadillo)

add_executable(predictSigmaGammaProton_SoftPomeron.exe physics/predictSigmaGammaProton_SoftPomeron.cpp)
target_link_libraries(predictSigmaGammaProton_SoftPomeron.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaProtonProton_HardPomeron.exe physics/fitSigmaProtonProton_HardPomeron.cpp)
target_link_libraries(fitSigmaProtonProton_HardPomeron.exe ihqcd_lib blas armadillo)

add_executable(fitSigmaProtonProton_SoftPomeron.exe physics/fitSigmaProtonProton_SoftPomeron.cpp)
target_link_libraries(fitSigmaProtonProton_SoftPomeron.exe ihqcd_lib blas armadillo)

add_executable(predictSigmaProtonProton_HardPomeron.exe physics/predictSigmaProtonProton_HardPomeron.cpp)
target_link_libraries(predictSigmaProtonProton_HardPomeron.exe ihqcd_lib blas armadillo)

add_executable(predictSigmaProtonProton_SoftPomeron.exe physics/predictSigmaProtonProton_SoftPomeron.cpp)
target_link_libraries(predictSigmaProtonProton_SoftPomeron.exe ihqcd_lib blas armadillo)

add_executable(globalFit_fixedKernel.exe physics/globalFit_fixedKernel.cpp)
target_link_libraries(globalFit_fixedKernel.exe ihqcd_lib blas armadillo)

add_executable(globalPomeronFit.exe physics/globalPomeronFit.cpp)
target_link_libraries(globalPomeronFit.exe ihqcd_lib blas armadillo)

add_executable(fitF2PhotonSigmaGammaGammaDISSigmaGammaProton.exe physics/fitF2PhotonSigmaGammaGammaDISSigmaGammaProton.cpp)
target_link_libraries(fitF2PhotonSigmaGammaGammaDISSigmaGammaProton.exe ihqcd_lib blas armadillo)

add_executable(gravitational_couplings_aux_calculations.exe physics/gravitational_couplings_aux_calculations.cpp)
target_link_libraries(gravitational_couplings_aux_calculations.exe ihqcd_lib blas armadillo)

add_executable(test_f2Spectrum.exe tests/test_f2Spectrum.cpp)
target_link_libraries(test_f2Spectrum.exe ihqcd_lib blas armadillo)

add_executable(fitGammaGammaAllPars.exe physics/fitGammaGammaAllPars.cpp)
target_link_libraries(fitGammaGammaAllPars.exe ihqcd_lib blas armadillo)

add_executable(fitGammaProtonAllPars.exe physics/fitGammaProtonAllPars.cpp)
target_link_libraries(fitGammaProtonAllPars.exe ihqcd_lib blas armadillo)

add_executable(fitGammaGammaGammaProtonAllPars.exe physics/fitGammaGammaGammaProtonAllPars.cpp)
target_link_libraries(fitGammaGammaGammaProtonAllPars.exe ihqcd_lib blas armadillo)

add_executable(fitGammaGammaKernelPars.exe physics/fitGammaGammaKernelPars.cpp)
target_link_libraries(fitGammaGammaKernelPars.exe ihqcd_lib blas armadillo)

add_executable(fitGammaProtonKernelPars.exe physics/fitGammaProtonKernelPars.cpp)
target_link_libraries(fitGammaProtonKernelPars.exe ihqcd_lib blas armadillo)

add_executable(fitGammaGammaGammaProtonKernelPars.exe physics/fitGammaGammaGammaProtonKernelPars.cpp)
target_link_libraries(fitGammaGammaGammaProtonKernelPars.exe ihqcd_lib blas armadillo)