#include <iostream>
#include <string>
#include <vector>
#include "F2Photon.h"
#include "HardPomeron.h"
#include "schrodinger/schrodinger.h"
#include "methods/search.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    string data_path = "expdata/F2_photon/F2_photon_xmax_0.01.txt";
    if (argc == 2) data_path = argv[1]; 

    F2Photon f2(data_path);

    // Testing if experimental data is dealt correctly
    vector<double> f2vals = f2.expVal(), errors = f2.expErr();
    vector<vector<double> > kinPts = f2.expKinematics();
    cout << "Q2\tx\tF2\terror" << endl;
    for(int i = 0; i < f2vals.size(); i++) cout << kinPts[0][i] << '\t' << kinPts[1][i] << '\t' << f2vals[i] << '\t' << errors[i] << endl;

    // Compute Chebyschev matrices
    chebSetN(800);

    // Setup gluon kernel and compute the Reggeons for t = 0
    HardPomeron hard;
    hard.computeReggeTrajectories();
    vector<Reggeon> reggeons = computeReggeons(hard, 0.0, 4);
    Spectra spec(0.0, reggeons);
    // Compute the IzNs using the IzN function
    cout << "Testing IzN function" << endl;
    cout << "Q2\tx\tIzN.1\tIzN.2\tIzN.3\tIzN.4" << endl;
    for(int i = 0; i < kinPts[0].size(); i++) 
    {
        cout << kinPts[0][i] << '\t' << kinPts[1][i] << '\t' << f2.IzN({kinPts[0][i]}, reggeons[0]) << '\t'; 
        cout << f2.IzN({kinPts[0][i]}, reggeons[1]) << '\t' << f2.IzN({kinPts[0][i]}, reggeons[2]) << '\t' ;
        cout << f2.IzN({kinPts[0][i]}, reggeons[3]) << endl;
    }
    // Compute the IzNs using the getIzs function
    cout << "Testing getIzs function" << endl;
    vector<kinStruct> izs = f2.getIzs(kinPts, {spec});
    cout << "Expecting izs to be of size 1" << endl;
    cout << "Size of izs: " << izs.size() << endl;
    vector<double> izns = izs[0].izns;
    cout << "Q2\tIzN.1\tIzN.2\tIzN.3\tIzN.4" << endl;
    for(int i = 0; i < izs.size(); i++)
    {
        cout << izs[i].kinVar << '\t';
        vector<double> izns = izs[i].izns;
        for(int j = 0; j < izns.size(); j++) cout << izns[j] << '\t';
        cout << endl;
    }
    // Compute the IzNBars using the IzNBar function
     // gs vector
    vector<double> gs = {0, 0, 0, 0};
    cout << "j0: " << reggeons[0].getJ() << " j1: " << reggeons[1].getJ() << " j2: " << reggeons[2].getJ() << " j3: " << reggeons[3].getJ() << endl;
    cout << "Testing IzNBar function" << endl;
    cout << "x\tIzNBar.1\tIzNBar.2\tIzNBar.3\tIzNBar.4" << endl;
    for(int i = 0; i < kinPts[1].size(); i++)
    {
        cout << kinPts[1][i] << '\t' << f2.IzNBar({kinPts[1][i]}, reggeons[0], gs) << '\t' << f2.IzNBar({kinPts[1][i]}, reggeons[1], gs) << '\t';
        cout << f2.IzNBar({kinPts[1][i]}, reggeons[2], gs) << '\t' << f2.IzNBar({kinPts[1][i]}, reggeons[3], gs) << endl;
    }
    // Compute the IzNBars using the function getIzsBar
    vector<kinStruct> izbars = f2.getIzsBar(kinPts, {spec}, gs);
    cout << "Testing getIzsBar function" << endl;
    cout << "Check that we have the same values" << endl;
    for(int i = 0; i < kinPts[1].size(); i++)
    {
        kinStruct iznbar = binary_search<kinStruct>(izbars, kinStruct(kinPts[1][i], {}));
        vector<double> iznbars = iznbar.izns;
        cout << kinPts[1][i] << '\t';
        for(int i = 0; i < iznbars.size(); i++) cout << iznbars[i] << '\t';
        cout << endl;
    }

    cout << "Testing predict function" << endl;
    vector<double> preds = f2.predict(izs, izbars, kinPts, false);

    double chi2 = f2.rss(izs, izbars, kinPts);
    cout << "chi2: " << chi2 << endl;

    return 0;
}