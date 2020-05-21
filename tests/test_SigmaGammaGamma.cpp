#include <iostream>
#include <string>
#include <vector>
#include "SigmaGammaGamma.h"
#include "SoftPomeron.h"
#include "schrodinger/schrodinger.h"
#include "methods/search.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    string data_path = "expdata/gammagammaScatteringL3Processed.txt";
    if (argc == 2) data_path = argv[1];

    cout << "Loading data from " << data_path << endl; 

    SigmaGammaGamma sigma(data_path);

    // Testing if experimental data is dealt correctly
    vector<double> sigmas = sigma.expVal(), sigmaErr = sigma.expErr();
    vector<vector<double> > Ws = sigma.expKinematics();
    cout << "W(GeV)\tWPlus\tWMinus\tsigma\tsigmaErr" << endl;
    for(int i = 0; i < sigmas.size(); i++) cout << Ws[0][i] << '\t' << Ws[1][i] << '\t' << Ws[2][i] << '\t' << sigmas[i] << '\t' << sigmaErr[i] << endl;

    // Compute Chebyschev matrices
    chebSetN(1000);

    // Setup gluon kernel and compute the Reggeons for t = 0
    SoftPomeron soft;
    soft.computeReggeTrajectories();
    vector<Reggeon> reggeons = computeReggeons(soft, 0.0, 2);
    Spectra spec(0.0, reggeons);

    // Compute the IzNs using the IzN function
    cout << "Testing IzN function" << endl;
    cout << "W\tIzN.1\tIzN.2" << endl;
    for(int i = 0; i < Ws[0].size(); i++) cout << Ws[0][i] << '\t' << sigma.IzN({Ws[0][i]}, reggeons[0]) << '\t' << sigma.IzN({Ws[0][i]}, reggeons[1]) << endl;
    // Compute the IzNs using the getIzs function
    cout << "Testing getIzs function" << endl;
    vector<kinStruct> izs = sigma.getIzs(Ws, {spec});
    cout << "Expecting izs to be of size 1" << endl;
    cout << "Size of izs: " << izs.size() << endl;
    vector<double> izns = izs[0].izns;
    cout << "IzN.1\tIzN.2" << endl;
    for(int i = 0; i < izns.size(); i++) cout << izns[i] << '\t';
    cout << endl;

    // Compute the IzNBars using the getIzsBar function
    // gs vector
    vector<double> gs = {1, 2};
    cout << "j0: " << reggeons[0].getJ() << " j1: " << reggeons[1].getJ() << endl;
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