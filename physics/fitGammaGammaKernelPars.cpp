#include <iostream>
#include <string>
#include <complex>
#include <random>
#include "IHQCD.h"
#include "F2Photon.h"
#include "SigmaGammaGamma.h"
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
    string f2photon_path, sigma_gg_path;
    if(argc < 12)
    {
        f2photon_path = "expdata/F2_photon/F2Photon_data.txt";
        sigma_gg_path = "expdata/GammaGamma/SigmaGammaGamma_PDG_data_W_gt_10.txt";
        Im_g0 = -0.000175226; Im_g1 = 0.000230161; Im_g2 = -0.000154943; Im_g3 = 0.0014227;
        invls = 6.491; a = -4.567; b = 1.485; c = 0.653; d = -0.113;
        cout << "Program usage: " + string(argv[0]) + " f2photon_path sigma_gg_path invls a b c d k1 k2 k3 k4" << endl;
        cout << "Using default values." << endl;
    }
    else
    {
        f2photon_path = argv[1]; sigma_gg_path = argv[2];
        invls = stod(argv[3]); a = stod(argv[4]); b = stod(argv[5]); c = stod(argv[6]); d = stod(argv[7]);
        Im_g0 = stod(argv[8]); Im_g1 = stod(argv[9]); Im_g2 = stod(argv[10]); Im_g3 = stod(argv[11]);
    }
    cout << "Starting fit with values:" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
    cout << "Using the following Im_gns values:" << endl;
    cout << "Im_g0: " << Im_g0 << " Im_g1: " << Im_g1 << " Im_g2: " << Im_g2 << " Im_g3: " << Im_g3 << endl;

    // Define the process observables and load the data needed for the fit
    F2Photon f2photon(f2photon_path);
    SigmaGammaGamma sigma_gg(sigma_gg_path);

    // Get the experimental points
    int npoints = 0;
    vector<vector<double> > f2photon_pts = f2photon.expKinematics();
    npoints += f2photon_pts[0].size();
    vector<vector<double> > sigma_gg_pts = sigma_gg.expKinematics();
    npoints += sigma_gg_pts[0].size();

    // Setup Chebyschev computation
    chebSetN(1000);

    // Setup HardPomeron Kernel and compute the Reggeons
    HardPomeron hard(4 , {invls, a, b, c, d});
    const vector<double> Im_gns = {Im_g0, Im_g1, Im_g2, Im_g3};

    // Definition of the functions we want to fit
    auto f = [&f2photon, &sigma_gg, &f2photon_pts, &sigma_gg_pts, &hard, &Im_gns] (const std::vector<double> & X)
    {
        // X is a vector with the values of invls, a, b, c, d, Im_g0, Im_g1, Im_g2, Im_g3
        vector<double> kernelPars = {X[0], X[1], X[2], X[3], X[4]};
        
        cout << "invls: " << X[0] << " a: " << X[1] << " b: " << X[2] << " c: " << X[3] << " d: " << X[4] << endl;
        
        hard.computeReggeTrajectories(kernelPars);
        vector<Reggeon> reggeons = computeReggeons(hard, 0, 4);
        Spectra spec(0, reggeons);
        vector<Spectra> spectrum = {spec};

        // Compute IzNs
        vector<kinStruct> F2Photon_IzNs = f2photon.getIzs(f2photon_pts, spectrum);
        vector<kinStruct> sigma_gg_IzNs = sigma_gg.getIzs(sigma_gg_pts, spectrum);

        // Compute the IzNsBar
        vector<kinStruct> F2Photon_IzNBars = f2photon.getIzsBar(f2photon_pts, spectrum, Im_gns);
        vector<kinStruct> sigma_gg_IzNBars = sigma_gg.getIzsBar(sigma_gg_pts, spectrum, Im_gns);

        // Compute the chi2 of each process
        double F2Photon_chi2 = f2photon.chi2(F2Photon_IzNs, F2Photon_IzNBars, f2photon_pts);
        double sigma_gg_chi2 = sigma_gg.chi2(sigma_gg_IzNs, sigma_gg_IzNBars, sigma_gg_pts);
        double chi2 = F2Photon_chi2 + sigma_gg_chi2;
        cout << "chi2: " << chi2 << endl; 
        return chi2;
    };


    // Start the fit now
    vector<double> X_guess = {invls, a, b, c, d};
    vector<double> deltas = {0.5, 0.5, 0.5, 0.5, 0.5};
    vector<double> X_opt = optimFunction(X_guess, f, deltas);
    
    // Print X_opt
    cout << "Found a minimum for: ";
    for(int i = 0; i < X_opt.size(); i++) cout << X_opt[i] << " ";
    cout << endl;
    cout << "Number of experimental points: " << npoints << endl;
    cout << "Number of degrees of freedom: " << npoints - X_opt.size() << endl;
    cout << "Best chi2 / Ndof: " << f(X_opt) / (npoints - X_opt.size()) << endl;


    return 0;
}