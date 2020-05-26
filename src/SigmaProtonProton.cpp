#include <iostream>
#include "SigmaProtonProton.h"

const double mb_to_GEVMINUS2 = 10 / 3.894 ;

SigmaProtonProton::SigmaProtonProton(std::string data_path): Sigma(data_path, mb_to_GEVMINUS2) 
{
    std::cout << "Loaded sigma(p p -> X) data." << std::endl;
}

double SigmaProtonProton::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    return 1;
}