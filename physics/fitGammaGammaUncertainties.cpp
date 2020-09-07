#include <iostream>
#include <vector>
#include <armadillo>
#include "IHQCD.h"
#include "F2Photon.h"
#include "SigmaGammaGamma.h"
#include "HardPomeron.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"
#include "uncertainties.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    arma::mat hessian;
    if(argc < 2)
    {
        cout << "Program usage: " << argv[0] << " f2_photon_data sigma_gg_data Im_g0 Im_g1 Im_g2 Im_g3 h" << endl;
        cout << "Or " << argv[0] << " hessian_path" << endl;
        return 0;
    }
    else if (argc == 2)
    {
        hessian.load(string(argv[1])); 
    }
    else if (argc == 8)
    {
        string data_path_f2 = argv[1], data_path_sigma = argv[2];
        double Im_g0 = stod(argv[3]), Im_g1 = stod(argv[4]), Im_g2 = stod(argv[5]), Im_g3 = stod(argv[6]);
        double h = stod(argv[7]);
        
        cout << "Computing uncertainties for values :" << endl;
        cout << "Im_g0: " << Im_g0 << " Im_g1: " << Im_g1 <<  " Im_g2: " << Im_g2 << " Im_g3: " << Im_g3 << endl;  
        // Setup sigma(gamma gamma -> hadrons) object
        F2Photon photon(data_path_f2);
        SigmaGammaGamma sigma(data_path_sigma);

        // Setup HardPomeron Kernel and GNs vector
        HardPomeron hard(4, {6.46892, -4.69919, 1.12825, 0.664399, -0.0982592});
        vector<double> GNs = {Im_g0, Im_g1, Im_g2, Im_g3};

        // Setup HQCDP object
        HQCDP hqcdp;
        hqcdp.addProcessObservable(photon);
        hqcdp.addProcessObservable(sigma);
        hqcdp.addKernel(hard);
        hqcdp.setGNs(GNs);
    
        // Compute the spectrum to check Reggeon properties
        chebSetN(1000);
        hqcdp.computeSpectrum();

        // We now compute the hessian matrix
        auto f = [&hqcdp] (const vector<double> & x)
        {
            // Update gns
            hqcdp.setGNs(x);
            // return chi2
            return hqcdp.chi2();
        };
        hessian = computeHessian(f, GNs, h);
        saveMatrix(hessian, "gamma_gamma_hessian.txt");
    }
    else
    {
        cout << "Invalid number of arguments. Call this program without arguments for instructions." << endl;
        return 0;
    }
    // We can now compute the parameter uncertainties
    vector<double> parUncert = parameterUncertainties(hessian);

    // Print the parameter values and the corresponding uncertainties
    cout << "Delta_Im_g0: " << parUncert[0] << endl;
    cout << "Delta_Im_g1: " << parUncert[1] << endl;
    cout << "Delta_Im_g2: " << parUncert[2] << endl;
    cout << "Delta_Im_g3: " << parUncert[3] << endl;

    return 0;

}