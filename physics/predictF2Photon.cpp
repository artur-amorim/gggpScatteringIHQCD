#include <iostream>
#include <fstream>
#include <string>
#include "HQCDP.h"
#include "IHQCD.h"
#include "HardPomeron.h"
#include "Reggeon.h"
#include "F2Photon.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2, g3, g4;
    string data_path_f2, data_path_sigma;

    if (argc < 6)
    {
        g1 = -0.000178883; g2 = 0.000193696; g3 = -0.000319052; g4 = 0.00133664;
        data_path_f2 = "expdata/F2_photon/F2Photon_data.txt";
    }
    else
    {
        data_path_f2 = argv[1];
        g1 = stod(argv[2]); g2 = stod(argv[3]); g3 = stod(argv[4]); g4 = stod(argv[5]);
    }

    cout << "Using values:" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;
    
    // Create F2 object
    F2Photon f2(data_path_f2);
    vector<vector<double> > F2points = f2.expKinematics();

    /// Setup HardPomeron Kernel and GNs vector
    HardPomeron hard(4, {6.46892, -4.69919, 1.12825, 0.664399, -0.0982592});
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addKernel(hard);
    hqcdp.setGNs(GNs);

    // Compute the spectrum to geth the reggeons
    chebSetN(1000);
    hqcdp.computeSpectrum();
    vector<Spectra> spectrum = hqcdp.getSpectrum();
    vector<Reggeon> reggeons = spectrum[0].getReggeons();

    // Compute IzNBars
    cout << "Computing F2Photon IzNBars" << endl;
    vector<kinStruct> F2IzNBars = f2.getIzsBar(F2points, spectrum, GNs);

    // Compute IzNs
    cout << "Computing F2Photon IzNs" << endl;
    vector<kinStruct> F2IzNs = f2.getIzs(F2points, spectrum);

    
    // Compute F2Photon chi2
    double F2chi2 = f2.chi2(F2IzNs, F2IzNBars, F2points);
    const int nPoints = F2points[0].size();
    cout << "Number of points: " << nPoints << endl;
    cout << "The F2 chi2 / points is " << F2chi2 / nPoints << endl;
    cout << "The F2 chi2 / Ndof is " << F2chi2 / (nPoints - 4) << endl;


    // Compute more predicted points in order to plot
    // Getting list of Q2s of F2
    std::vector<double> F2Q2s = f2.getDataPts()[0];
    sort(F2Q2s.begin(), F2Q2s.end() );
    F2Q2s.erase( unique( F2Q2s.begin(), F2Q2s.end() ), F2Q2s.end() );

    // Now we can generate the points where we want to compute the F2predictions
    // and also compute the predictions
    // Clear the necessary vectors
    F2points.clear(); F2IzNBars.clear();
    F2IzNs.clear();

    // Compute F2 points
    vector<double> Q2s(0), xs(0);
    for(int i = 0; i < F2Q2s.size(); i++)
    {
        for(double x = 1e-3; x <= 3e-2; x += 1e-6)
        {
            Q2s.push_back(F2Q2s[i]);
            xs.push_back(x);
        }
    }
    F2points = {Q2s, xs};
    
    cout << "Predicting F2Photon" << endl;

    // Predict the central values of F2
    F2IzNBars = f2.getIzsBar(F2points, spectrum, GNs);
    F2IzNs = f2.getIzs(F2points, spectrum);
    vector<double> F2pred = f2.predict(F2IzNs, F2IzNBars, F2points, true);

    return 0;
}