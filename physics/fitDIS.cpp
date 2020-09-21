#include <iostream>
#include <string>
#include "IHQCD.h"
#include "F2.h"
#include "FL.h"
#include "HardPomeron.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    // Gluon Kernel parameters
    double invls, a, b, c, d;
    // Gravitational photon couplings
    double Im_g0, Im_g1, Im_g2, Im_g3;
    string f2_path, fl_path;
    if(argc < 12)
    {
        f2_path = "expdata/DIS/F2_data.txt";
        fl_path = "expdata/DIS/FL_data.txt";
        Im_g0 = -0.0498758; Im_g1 = 0.0199297; Im_g2 = -0.00480236; Im_g3 = 0.349284;
        invls = 6.491; a = -4.567; b = 1.485; c = 0.653; d = -0.113;
        cout << "Program usage: " + string(argv[0]) + " f2_path fl_path invls a b c d k1 k2 k3 k4" << endl;
        cout << "Using default values." << endl;
    }
    else
    {
        f2_path = argv[1]; fl_path = argv[2];
        invls = stod(argv[3]); a = stod(argv[4]); b = stod(argv[5]); c = stod(argv[6]); d = stod(argv[7]);
        Im_g0 = stod(argv[8]); Im_g1 = stod(argv[9]); Im_g2 = stod(argv[10]); Im_g3 = stod(argv[11]);
    }
    cout << "Starting fit with values:" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
    cout << "Im_g0: " << Im_g0 << " Im_g1: " << Im_g1 << " Im_g2: " << Im_g2 << " Im_g3: " << Im_g3 << endl;

    // Define the process observables and load the data needed for the fit
    F2 f2(f2_path);
    FL fl(fl_path);

    // Get the experimental points
    int npoints = 0;
    vector<vector<double> > f2_pts = f2.expKinematics();
    npoints += f2_pts[0].size();
    vector<vector<double> > fl_pts = fl.expKinematics();
    npoints += fl_pts[0].size();

    // Setup Chebyschev computation
    chebSetN(1000);

    // Setup HardPomeron Kernel and compute the Reggeons
    HardPomeron hard(4 , {invls, a, b, c, d});


    // Definition of the functions we want to fit
    auto f = [&f2, &fl, &f2_pts, &fl_pts, &hard] (const std::vector<double> & X)
    {
        // X is a vector with the values of invls, a, b, c, d, Im_g0, Im_g1, Im_g2, Im_g3
        vector<double> kernelPars = {X[0], X[1], X[2], X[3], X[4]};
        vector<double> Im_gns = {X[5], X[6], X[7], X[8]};

        cout << "invls: " << X[0] << " a: " << X[1] << " b: " << X[2] << " c: " << X[3] << " d: " << X[4] << endl;
        cout << "Im_g0: " << X[5] << " Im_g1: " << X[6] << " Im_g2: " << X[7] << " Im_g3: " << X[8] << endl;

        hard.computeReggeTrajectories(kernelPars);
        vector<Reggeon> reggeons = computeReggeons(hard, 0, 4);
        Spectra spec(0, reggeons);
        vector<Spectra> spectrum = {spec};

        // Compute IzNs
        vector<kinStruct> F2_IzNs = f2.getIzs(f2_pts, spectrum);
        vector<kinStruct> FL_IzNs = fl.getIzs(fl_pts, spectrum);

        // Compute the IzNsBar
        vector<kinStruct> F2_IzNBars = f2.getIzsBar(f2_pts, spectrum, Im_gns);
        vector<kinStruct> FL_IzNBars = fl.getIzsBar(fl_pts, spectrum, Im_gns);

        // Compute the chi2 of each process
        double F2_chi2 = f2.chi2(F2_IzNs, F2_IzNBars, f2_pts);
        double FL_chi2 = fl.chi2(FL_IzNs, FL_IzNBars, fl_pts);
        double chi2 = F2_chi2 + FL_chi2;
        cout << "chi2: " << chi2 << endl; 
        return chi2;
    };


    // Start the fit now
    vector<double> X_guess = {invls, a, b, c, d, Im_g0, Im_g1, Im_g2, Im_g3};
    vector<double> deltas = {0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1};
    vector<double> X_opt = optimFunction(X_guess, f, deltas);
    
    // Print X_opt
    cout << "Found a minimum for: ";
    for(int i = 0; i < X_opt.size(); i++) cout << X_opt[i] << " ";
    cout << endl;
    cout << "Number of experimental points: " << npoints << endl;
    cout << "Number of degrees of freedom: " << npoints - X_opt.size() << endl;
    double chi2 = f(X_opt);
    cout << "Best chi2: " << chi2 << endl;
    cout << "Best chi2 / Ndof: " << chi2 / (npoints - X_opt.size()) << endl;


    return 0;
}