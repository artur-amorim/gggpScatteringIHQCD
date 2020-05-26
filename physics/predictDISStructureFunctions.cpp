#include <iostream>
#include <fstream>
#include <string>
#include "HQCDP.h"
#include "IHQCD.h"
#include "HardPomeron.h"
#include "Reggeon.h"
#include "F2.h"
#include "FL.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double g1, g2, g3, g4;

    if (argc < 5)
    {
        g1 = 0; g2 = 0; g3 = 0; g4 = 0;
        cout << "Using default values:" << endl;
    }
    else
    {
        g1 = stod(argv[1]); g2 = stod(argv[2]); g3 = stod(argv[3]); g4 = stod(argv[4]);
        cout << "Using values:" << endl;
    }

    cout << " g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl; 
    
    // Create F2 object
    F2 f2("expdata/DIS/F2_data.txt");
    vector<vector<double> > F2points = f2.expKinematics();
    FL fl("expdata/DIS/FL_data.txt");
    vector<vector<double> > FLpoints = fl.expKinematics();

    // Setup HQCDP kernel and compute the reggeons
    HQCDP hqcdp;
    // Add kernel
    HardPomeron hard;
    hqcdp.addKernel(hard);
    // Add f2 and fl to the list of processes of HQCDP
    hqcdp.addProcessObservable(f2);
    hqcdp.addProcessObservable(fl);
    // Add the gns
    vector<double> GNs = {g1, g2, g3, g4};
    hqcdp.setGNs(GNs);

    // Compute the spectrum and get the reggeons
    chebSetN(1000);
    hqcdp.computeSpectrum();
    vector<Spectra> spectrum = hqcdp.getSpectrum();
    vector<Reggeon> reggeons = spectrum[0].getReggeons();

    // Compute IzNBars
    cout << "Computing F2 IzNBars" << endl;
    vector<kinStruct> F2IzNBars = f2.getIzsBar(F2points, spectrum, GNs);
    cout << "Computing FL IzNBars" << endl;
    vector<kinStruct> FLIzNBars = fl.getIzsBar(FLpoints, spectrum, GNs);

    // Compute IzNs
    cout << "Computing F2 IzNs" << endl;
    vector<kinStruct> F2IzNs = f2.getIzs(F2points, spectrum);
    cout << "Computing FL IzNs" << endl;
    vector<kinStruct> FLIzNs = fl.getIzs(FLpoints, spectrum);

    // Compute F2 and FL
    std::cout << "Predicting F2 for the given values of Q2 and x" << std::endl;
    vector<double> F2pred = f2.predict(F2IzNs, F2IzNBars, F2points, false);
    std::cout << "Predicting FL for the given values of Q2 and x" << std::endl;
    vector<double> FLpred = fl.predict(FLIzNs, FLIzNBars, FLpoints, false);
    
    // Compute F2 and FL chi2
    double F2chi2 = f2.chi2(F2IzNs, F2IzNBars, F2points);
    cout << "The F2 chi2 is " << F2chi2 / (F2points[0].size()) << endl;
    double FLchi2 = fl.chi2(FLIzNs, FLIzNBars, FLpoints);
    cout << "The FL chi2 is " << FLchi2 / FLpoints[0].size() << endl;
    cout << "Total chi2 is " << (F2chi2 + FLchi2) / (F2points[0].size() + FLpoints[0].size()) << endl;
    // Compute more predicted points in order to plot
    // Getting list of Q2s of F2
    std::vector<double> F2Q2s = f2.getDataPts()[0];
    sort(F2Q2s.begin(), F2Q2s.end() );
    F2Q2s.erase( unique( F2Q2s.begin(), F2Q2s.end() ), F2Q2s.end() );
    // Getting list of Q2s of FL
    std::vector<double> FLQ2s = fl.getDataPts()[0];
    sort(FLQ2s.begin(), FLQ2s.end() );
    FLQ2s.erase( unique( FLQ2s.begin(), FLQ2s.end() ), FLQ2s.end() );
    // Now we can generate the points where we want to compute the F2 and FL predictions
    // and also compute the predictions
    // Clear the necessary vectors
    F2points.clear(); FLpoints.clear();
    F2IzNBars.clear(); FLIzNBars.clear();
    F2IzNs.clear(); FLIzNs.clear();
    F2pred.clear(); FLpred.clear();

    // Compute F2 points
    vector<double> Q2s(0), xs(0);
    for(int i = 0; i < F2Q2s.size(); i++)
    {
        for(double x = 1e-6; x <= 1e-2; x += 1e-6)
        {
            Q2s.push_back(F2Q2s[i]);
            xs.push_back(x);
        }
    }
    F2points = {Q2s, xs};
    // Compute FL points
    Q2s.clear(); xs.clear();
    for(int i = 0; i < FLQ2s.size(); i++)
    {
        for(double x = 1e-6; x <= 1e-2; x += 1e-6)
        {
            Q2s.push_back(FLQ2s[i]);
            xs.push_back(x);
        }
    }
    FLpoints = {Q2s, xs};
    
    cout << "Predicting DIS Structure Functions F2 and FL." << endl;

    // Predict the central values of F2
    F2IzNBars = f2.getIzsBar(F2points, spectrum, GNs);
    F2IzNs = f2.getIzs(F2points, spectrum);
    F2pred = f2.predict(F2IzNs, F2IzNBars, F2points, true);
    // Predict the central values of FL
    FLIzNBars = fl.getIzsBar(FLpoints, spectrum, GNs);
    FLIzNs = fl.getIzs(FLpoints, spectrum);
    FLpred = fl.predict(FLIzNs, FLIzNBars, FLpoints, true);
 
    return 0;
}