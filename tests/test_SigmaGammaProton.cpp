#include <iostream>
#include <string>
#include <vector>
#include "SigmaGammaProton.h"
#include "HardPomeron.h"
#include "schrodinger/schrodinger.h"
#include "methods/search.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    string data_path = "expdata/GammaP/SigmaGammaProton.txt";
    double g1, g2, g3, g4;
    if (argc < 6)
    {
        data_path = "expdata/GammaP/SigmaGammaProton.txt";
        g1 = 0.0; g2 = 0.0; g3 = 0.0; g4 = 0.0;
    }
    else
    {
        data_path = argv[1];
        g1 = stod(argv[2]); g2 = stod(argv[3]); g3 = stod(argv[4]); g4 = stod(argv[5]);
    }
    

    cout << "Loading data from " << data_path << endl; 

    SigmaGammaProton sigma(data_path);

    // Testing if experimental data is dealt correctly
    vector<double> sigmas = sigma.expVal(), sigmaErr = sigma.expErr();
    vector<vector<double> > Ws = sigma.expKinematics();
    cout << "W(GeV)\tWPlus\tWMinus\tsigma\tsigmaErr" << endl;
    for(int i = 0; i < sigmas.size(); i++) cout << Ws[0][i] << '\t' << Ws[1][i] << '\t' << Ws[2][i] << '\t' << sigmas[i] << '\t' << sigmaErr[i] << endl;

    // Compute Chebyschev matrices
    chebSetN(1000);

    // Setup gluon kernel and compute the Reggeons for t = 0
    HardPomeron hard;
    hard.computeReggeTrajectories();
    vector<Reggeon> reggeons = computeReggeons(hard, 0.0, 4);
    Spectra spec(0.0, reggeons);

    // Compute the IzNs using the IzN function
    cout << "Testing IzN function" << endl;
    cout << "W\tIzN.1\tIzN.2\tIzN.3\tIzN.4" << endl;
    for(int i = 0; i < Ws[0].size(); i++)
    {
        cout << Ws[0][i] << '\t';
        for(int j = 0; j < reggeons.size(); j++) cout << sigma.IzN({Ws[0][i]}, reggeons[j]) << '\t';
        cout << endl;
    }
    // Compute the IzNs using the getIzs function
    cout << "Testing getIzs function" << endl;
    vector<kinStruct> izs = sigma.getIzs(Ws, {spec});
    cout << "Expecting izs to be of size 1" << endl;
    cout << "Size of izs: " << izs.size() << endl;
    vector<double> izns = izs[0].izns;
    cout << "IzN.1\tIzN.2\tIzN.3\tIzN.4" << endl;
    for(int i = 0; i < izns.size(); i++) cout << izns[i] << '\t';
    cout << endl;

    // Compute the IzNBars using the getIzsBar function
    // gs vector
    vector<double> gs = {g1, g2, g3, g4};
    for(int i = 0; i < reggeons.size(); i++) cout << "j"+to_string(i)+": " << reggeons[i].getJ() << '\t';
    cout << endl;
    // Compute the IzNBars using the function getIzsBar
    vector<kinStruct> izbars = sigma.getIzsBar(Ws, {spec}, gs);
    cout << "Testing getIzsBar function" << endl;
    cout << "izbars should have the same size as Ws[0]" << endl;
    cout << "izbars size: " << izbars.size() << " Ws[0] size: " << Ws[0].size() << endl;
    for(int i = 0; i < Ws[0].size(); i++)
    {
        kinStruct iznbar = binary_search<kinStruct>(izbars, kinStruct(Ws[0][i], {}));
        vector<double> iznbars = iznbar.izns;
        cout << Ws[0][i] << '\t';
        for(int i = 0; i < iznbars.size(); i++) cout << iznbars[i] << '\t';
        cout << endl;
    }

    cout << "Testing predict function" << endl;
    vector<double> preds = sigma.predict(izs, izbars, Ws, true);
    
    // Now let's test diffObsWeighted
    vector<double> obs_weighted = sigma.diffObsWeighted(izs, izbars, Ws);

    double chi2 = sigma.chi2(izs, izbars, Ws);
    cout << "chi2: " << chi2 << endl;

    return 0;
}