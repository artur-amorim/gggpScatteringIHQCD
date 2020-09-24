#include <iostream>
#include <fstream>
#include <string>
#include "HQCDP.h"
#include "IHQCD.h"
#include "HardPomeron.h"
#include "Reggeon.h"
#include "SigmaGammaProton.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double invls, a, b, c, d;
    double g1, g2, g3, g4;
    string data_path;
    if (argc < 11)
    {
        data_path = "expdata/GammaP/SigmaGammaP_PDG_data_W_gt_461.txt";
        invls = 6.46892; a = -4.69919; b = 1.12825; c = 0.664399; d = -0.0982592;
        g1 = -0.0510176; g2 = 0.017369; g3 = -0.0744977; g4 = 0.357739;
    }
    else
    {
        data_path = argv[1];
        invls = stod(argv[2]); a = stod(argv[3]); b = stod(argv[4]); c = stod(argv[5]); d = stod(argv[6]);
        g1 = stod(argv[7]); g2 = stod(argv[8]); g3 = stod(argv[9]); g4 = stod(argv[10]);
    }

    cout << "invls : " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
    cout << " g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl; 
    
    // Create SigmaGammaProton object
    SigmaGammaProton sigma(data_path);
    vector<vector<double> > sigma_points = sigma.expKinematics();

    /// Setup HardPomeron Kernel and GNs vector
    HardPomeron hard(4, {invls, a, b, c, d});
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(hard);
    hqcdp.setGNs(GNs);

    // Compute the spectrum to geth the reggeons
    chebSetN(1000);
    hqcdp.computeSpectrum();
    vector<Spectra> spectrum = hqcdp.getSpectrum();

    // Compute IzNBars
    cout << "Computing sigma(gamma p -> hadrons) IzNBars" << endl;
    vector<kinStruct> SigmaGammaPIzNBars = sigma.getIzsBar(sigma_points, spectrum, GNs);

    // Compute IzNs
    cout << "Computing sigma(gamma p -> hadrons) IzNs" << endl;
    vector<kinStruct> SigmaGammaPIzNs = sigma.getIzs(sigma_points, spectrum);
    
    // Compute sigma(gamma p -> hadrons) chi2
    const int nPoints = sigma_points[0].size();
    double sigma_chi2 = sigma.chi2(SigmaGammaPIzNs, SigmaGammaPIzNBars, sigma_points);
    cout << "Number of points: " << nPoints << endl;
    cout << "sigma(gamma p -> hadrons) chi2: " << sigma_chi2 << endl;
    cout << "The sigma(gamma p -> hadrons) chi2 / points is " << sigma_chi2 / nPoints << endl;
    cout << "The sigma(gamma p -> hadrons) chi2 / Ndof is " << sigma_chi2 / (nPoints -4 ) << endl;

    vector<double> Ws;
    for(double W = 1.5; W < 300; W += 0.1) Ws.push_back(W);
    vector<double> WPlus(Ws.size(), 0.0), WMinus(Ws.size(), 0.0);
    vector<vector<double>> kinPts = {Ws, WPlus, WMinus};
    // Compute new IzNBars
    SigmaGammaPIzNBars = sigma.getIzsBar(kinPts, spectrum, GNs);
    // Compute new IzNs
    SigmaGammaPIzNs = sigma.getIzs(kinPts, spectrum);
    // Compute  predictions of sigma(gamma p -> hadrons)
    std::cout << "Predicting sigma(gamma p -> hadrons) for the given values of sqrt(s)" << std::endl;
    vector<double> sigma_pred = sigma.predict(SigmaGammaPIzNs, SigmaGammaPIzNBars, kinPts, true);

    return 0;
}