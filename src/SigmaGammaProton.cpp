#include "SigmaGammaProton.h"

const double mub_to_GEVMINUS2 = 1.0 / (3.894e2) ;

SigmaGammaProton::SigmaGammaProton(std::string data_path): Sigma(data_path, mub_to_GEVMINUS2) {}