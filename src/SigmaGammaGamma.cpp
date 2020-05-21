#include "SigmaGammaGamma.h"

const double nb_to_GEVMINUS2 = 1.0 / (3.894e5) ;

SigmaGammaGamma::SigmaGammaGamma(std::string data_path): Sigma(data_path, nb_to_GEVMINUS2) {}