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
        g1 = 2.6527e-05; g2 = 0.000229015; g3 = 0.000155743; g4 = -0.00101869;
        data_path_f2 = "expdata/F2_photon/F2_photon_xmax_0.01.txt";
    }
    else
    {
        g1 = stod(argv[1]); g2 = stod(argv[2]); g3 = stod(argv[3]); g4 = stod(argv[4]);
        data_path_f2 = argv[5];
    }

    cout << "Using values:" << endl;
        cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;
    
    // Create F2 object
    F2Photon f2(data_path_f2);
    vector<vector<double> > F2points = f2.expKinematics();

    /// Setup HardPomeron Kernel and GNs vector
    HardPomeron hard;
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addKernel(hard);
    hqcdp.setGNs(GNs);

    // Compute the spectrum to geth the reggeons
    chebSetN(800);
    hqcdp.computeSpectrum();
    vector<Spectra> spectrum = hqcdp.getSpectrum();
    vector<Reggeon> reggeons = spectrum[0].getReggeons();

    // Compute IzNBars
    cout << "Computing F2 IzNBars" << endl;
    vector<kinStruct> F2IzNBars = f2.getIzsBar(F2points, spectrum, GNs);

    // Compute IzNs
    cout << "Computing F2 IzNs" << endl;
    vector<kinStruct> F2IzNs = f2.getIzs(F2points, spectrum);

    // Compute F2 and FL
    std::cout << "Predicting F2 for the given values of Q2 and x" << std::endl;
    vector<double> F2pred = f2.predict(F2IzNs, F2IzNBars, F2points, false);
    
    // Compute F2 and FL chi2
    double F2chi2 = f2.rss(F2IzNs, F2IzNBars, F2points);
    cout << "The F2 chi2 is " << F2chi2 / (F2points[0].size() -4 ) << endl;

    // Compute more predicted points in order to plot
    // Getting list of Q2s of F2
    std::vector<double> F2Q2s = f2.getDataPts()[0];
    sort(F2Q2s.begin(), F2Q2s.end() );
    F2Q2s.erase( unique( F2Q2s.begin(), F2Q2s.end() ), F2Q2s.end() );

    // Now we can generate the points where we want to compute the F2 and FL predictions
    // and also compute the predictions
    // Clear the necessary vectors
    F2points.clear(); F2IzNBars.clear();
    F2IzNs.clear(); F2pred.clear();

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
    
    cout << "Predicting Photon F2" << endl;

    // Predict the central values of F2
    F2IzNBars = f2.getIzsBar(F2points, spectrum, GNs);
    F2IzNs = f2.getIzs(F2points, spectrum);
    F2pred = f2.predict(F2IzNs, F2IzNBars, F2points, false);

    // Ok, now i can write all of this in a file
    ofstream myfile;
    myfile.open("F2Photon_pred.txt");
    myfile << "Q2\tx\tF2" << endl;
    for(int i = 0; i < F2points[0].size(); i++)
        myfile << F2points[0][i] << '\t' << F2points[1][i] << '\t' << F2pred[i] << endl;
    myfile.close();
    return 0;
}