#include <iostream>
#include <vector>
#include "Reggeon.h"
#include "Spectra.h"
#include "PhotonScattering.h"

using namespace std;

int main()
{
    PhotonScattering photon;

    // Check getDataPts and setDataPts
    cout << "Number of data points should be zero" << endl;
    vector<vector<double> > data_pts = photon.getDataPts();
    cout << "Number of data points: " << data_pts.size() << endl;
    cout << "Setting data points to {{1,2}, {3,4}, {5,6}}" << endl;
    photon.setDataPts({{1,2}, {3,4}, {5,6}});
    data_pts = photon.getDataPts();
    for(int i = 0; i < data_pts[0].size(); i++) cout << data_pts[0][i] << '\t' << data_pts[1][i] << '\t' << data_pts[2][i] << endl;

    // Testing expKinematics, expValues and expErros. The lists shoud be empty
    vector<vector<double> > kinematics = photon.expKinematics();
    vector<double> values = photon.expVal(), errors = photon.expErr();
    cout << "Size of kinematics, values and errors vectors: " << kinematics.size() << " " << values.size() << " " << errors.size() << endl;

    // Check getNeededTVals function
    vector<double> tvals = photon.getNeededTVals();
    cout << "Testing getNeededTVals" << endl;
    cout << "Needed t values: ";
    for(int i = 0; i < tvals.size(); i++) cout << tvals[i] << '\t' << endl;

    // Testing IzN, IzNBar, getIzs and getIzsBar
    Reggeon reg;
    double izn = photon.IzN({}, reg), iznbar = photon.IzNBar({}, reg, {});
    cout << "IzN result: " << izn << " IzNBar result: " << iznbar << endl;
    vector<kinStruct> izs = photon.getIzs({{}}, {}), izbars = photon.getIzsBar({{}}, {}, {});
    cout << "Expected size of izs and izbars: 0" << endl;
    cout << "Size of izs: " << izs.size() << " Size of izbars: " << izbars.size() << endl;

    cout << "Expected size of predictions: 0" << endl;
    vector<double> predictions = photon.predict({}, {}, {{}}, false);
    cout << "predictions size: " << predictions.size() << endl;

    cout << "Expected size of diffObsWeigted: 0" << endl;
    vector<double> obs_weighted = photon.diffObsWeighted({}, {}, {{}});
    cout << "obs_weighted size: " << obs_weighted.size() << endl;

    cout << "Expected value of chi2: 0" << endl;
    double chi2 = photon.rss({}, {}, {{}});
    cout << "Chi2 obtained: " << chi2 << endl;

    // Testing assignment operator
    PhotonScattering photon2;
    photon2 = photon;

    cout << "If the assignment operator is working properly the output of getDataPts should be {{1, 2}, {3, 4}, {5, 6}}" << endl;
    data_pts = photon2.getDataPts();
    for(int i = 0; i < data_pts[0].size(); i++) cout << data_pts[0][i] << '\t' << data_pts[1][i] << '\t' << data_pts[2][i] << endl;

    // Testing copy constructor
    cout << "If the copy constructor is working properly the output of getDataPts should be {{1, 2}, {3, 4}, {5, 6}}" << endl;
    PhotonScattering photon3(photon);
    data_pts = photon3.getDataPts();
    for(int i = 0; i < data_pts[0].size(); i++) cout << data_pts[0][i] << '\t' << data_pts[1][i] << '\t' << data_pts[2][i] << endl;

    return 0 ;
}