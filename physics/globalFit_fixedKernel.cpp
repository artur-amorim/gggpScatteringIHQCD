#include <iostream>
#include <string>
#include <complex>
#include <random>
#include "IHQCD.h"
#include "F2Photon.h"
#include "F2.h"
#include "FL.h"
#include "SigmaGammaProton.h"
#include "SigmaGammaGamma.h"
#include "SigmaProtonProton.h"
#include "HardPomeron.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    // Gravitational photon couplings
    double k1, k2, k3, k4;
    // Gravitational proton couplings
    double kbar1, kbar2, kbar3, kbar4;
    string f2photon_path, f2_path, fl_path, sigma_gp_path, sigma_gg_path, sigma_pp_path;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    if(argc < 15)
    {
        f2photon_path = "expdata/F2_photon/F2Photon_data.txt";
        f2_path = "expdata/DIS/F2_data.txt";
        fl_path = "expdata/DIS/FL_data.txt";
        sigma_gp_path = "expdata/GammaP/SigmaGammaProton_2.txt";
        sigma_gg_path = "expdata/GammaGamma/gammagammaScatteringL3PHOJET.txt";
        sigma_pp_path = "expdata/SigmaProtonProton/SigmaProtonProton_data.txt";
        normal_distribution<> d_k1(0,1), d_k2(0, 1), d_k3(0, 1), d_k4(0, 1);
        normal_distribution<> d_kbar1(0,1), d_kbar2(0, 1), d_kbar3(0, 1), d_kbar4(0, 1);
        k1 = d_k1(gen); k2 = d_k2(gen); k3 = d_k3(gen); k4 = d_k4(gen);
        kbar1 = d_kbar1(gen); kbar2 = d_kbar2(gen); kbar3 = d_kbar3(gen); kbar4 = d_kbar4(gen);
        cout << "Program usage: " + string(argv[0]) + " f2photon_path f2_path fl_path sigma_gp_path sigma_gg_path sigma_pp_path k1 k2 k3 k4 ..." << endl;
        cout << "Using default values." << endl;
    }
    else
    {
        f2photon_path = argv[1]; f2_path = argv[2];
        fl_path = argv[3]; sigma_gp_path = argv[4];
        sigma_gg_path = argv[5]; sigma_pp_path = argv[6];
        normal_distribution<> d_k1(stod(argv[7]),0.1), d_k2(stod(argv[8]), 0.1), d_k3(stod(argv[9]), 0.1), d_k4(stod(argv[10]), 0.1);
        normal_distribution<> d_kbar1(stod(argv[11]), 0.1), d_kbar2(stod(argv[12]), 0.1), d_kbar3(stod(argv[13]), 0.1), d_kbar4(stod(argv[14]), 0.1);
        k1 = d_k1(gen); k2 = d_k2(gen); k3 = d_k3(gen); k4 = d_k4(gen);
        kbar1 = d_kbar1(gen); kbar2 = d_kbar2(gen); kbar3 = d_kbar3(gen); kbar4 = d_kbar4(gen);
    }
    cout << "Starting fit with values:" << endl;
    cout << "k1: " << k1 << " k2: " << k2 << " k3: " << k3 << " k4: " << k4 << endl;
    cout << "kbar1: " << kbar1 << " kbar2: " << kbar2 << " kbar3: " << kbar3 << " kbar4: " << kbar4 << endl;

    // Define the process observables and load the data needed for the fit
    F2Photon f2photon(f2photon_path);
    F2 f2(f2_path); FL fl(fl_path);
    SigmaGammaProton sigma_gp(sigma_gp_path);
    SigmaGammaGamma sigma_gg(sigma_gg_path);
    SigmaProtonProton sigma_pp(sigma_pp_path);

    // Get the experimental points
    int npoints = 0;
    vector<vector<double> > f2photon_pts = f2photon.expKinematics();
    npoints += f2photon_pts[0].size();
    vector<vector<double> > f2_pts = f2.expKinematics();
    npoints += f2_pts[0].size();
    vector<vector<double> > fl_pts = fl.expKinematics();
    npoints += fl_pts[0].size();
    vector<vector<double> > sigma_gp_pts = sigma_gp.expKinematics();
    npoints += sigma_gp_pts[0].size();
    vector<vector<double> > sigma_gg_pts = sigma_gg.expKinematics();
    npoints += sigma_gg_pts[0].size();
    vector<vector<double> > sigma_pp_pts = sigma_pp.expKinematics();
    npoints += sigma_pp_pts[0].size();

    // Setup Chebyschev computation
    chebSetN(1000);

    // Setup HardPomeron Kernel and compute the Reggeons
    HardPomeron hard;
    hard.computeReggeTrajectories();
    vector<Reggeon> reggeons = computeReggeons(hard, 0, 4);
    Spectra spec(0, reggeons);
    vector<Spectra> spectrum = {spec};

    // Ok now let's compute all the Izns, because our fitting parameters only enter in the IzNBars
    cout << "Computing F2Photon IzNs" << endl;
    vector<kinStruct> F2Photon_IzNs = f2photon.getIzs(f2photon_pts, spectrum);
    cout << "Computing F2 IzNs" << endl;
    vector<kinStruct> F2_IzNs = f2.getIzs(f2_pts, spectrum);
    cout << "Computing FL IzNs" << endl;
    vector<kinStruct> FL_IzNs = fl.getIzs(fl_pts, spectrum);
    cout << "Computing sigma(gamma p -> X) IzNs" << endl;
    vector<kinStruct> sigma_gp_IzNs = sigma_gp.getIzs(sigma_gp_pts, spectrum);
    cout << "Computing sigma(gamma gamma -> X) IzNs" << endl;
    vector<kinStruct> sigma_gg_IzNs = sigma_gg.getIzs(sigma_gg_pts, spectrum);
    cout << sigma_gg_IzNs.size() << '\t' << sigma_gg_IzNs[0].izns.size() << endl;
    cout << "Computing sigma(p p -> X) IzNs" << endl;
    vector<kinStruct> sigma_pp_IzNs = sigma_pp.getIzs(sigma_pp_pts, spectrum);

    auto f = [&f2photon, &f2, &fl, &sigma_gp, &sigma_gg, &sigma_pp, &spectrum,
             &f2photon_pts, &f2_pts, &fl_pts, &sigma_gp_pts, &sigma_gg_pts, &sigma_pp_pts,
             &F2Photon_IzNs, &F2_IzNs, &FL_IzNs, &sigma_gp_IzNs, &sigma_gg_IzNs, &sigma_pp_IzNs,
             &reggeons]
             (const std::vector<double> & X)
    {
        // X is a vector with the values of k1, k2, k3, k4, kbar1, kbar2, kbar3, kbar4
        vector<double> ks = {X[0], X[1], X[2], X[3]};
        vector<double> kbars = {X[4], X[5], X[6], X[7]};
        cout << "k1: " << ks[0] << " k2: " << ks[1] << " k3: " << ks[2] << " k4: " << ks[3] << " kbar1: " << kbars[0] << " kbar2: " << kbars[1] << " kbar3: " << kbars[2] << " kbar4: " << kbars[3] << endl;
        /*
            Now we need to compute the gn's according to the definition in the notes
        */
        vector<double> Im_gn_gg(reggeons.size()), Im_gn_gp(reggeons.size()), Im_gn_pp(reggeons.size());
        for(int i = 0; i < reggeons.size(); i++)
        {
            double jn = reggeons[i].getJ(), djndt = reggeons[i].getdJdt();
            // Common factor
            complex<double> gn(1/tan(M_PI_2 * jn), 1);
            gn = - M_PI_2 * gn * djndt / pow(2, jn);
            // Now we make the specific computations
            // The integrals that appear in the definition of Im_gn_gg are the same as the sigma_gg_IzNs
            Im_gn_gg[i] = imag(gn * ks[i]*ks[i] * sigma_gg_IzNs[0].izns[i]);
            Im_gn_gp[i] = imag(gn * ks[i] * kbars[i]);
            Im_gn_pp[i] = imag(gn * kbars[i]*kbars[i]);
        }

        // Compute the IzNsBar
        vector<kinStruct> F2Photon_IzNBars = f2photon.getIzsBar(f2photon_pts, spectrum, Im_gn_gg);
        vector<kinStruct> F2_IzNBars = f2.getIzsBar(f2_pts, spectrum, Im_gn_gp);
        vector<kinStruct> FL_IzNBars = fl.getIzsBar(fl_pts, spectrum, Im_gn_gp);
        vector<kinStruct> sigma_gp_IzNBars = sigma_gp.getIzsBar(sigma_gp_pts, spectrum, Im_gn_gp);
        vector<kinStruct> sigma_gg_IzNBars = sigma_gg.getIzsBar(sigma_gg_pts, spectrum, Im_gn_gg);
        vector<kinStruct> sigma_pp_IzNBars = sigma_pp.getIzsBar(sigma_pp_pts, spectrum, Im_gn_pp);

        // Compute the chi2 of each process
        double F2Photon_chi2 = f2photon.chi2(F2Photon_IzNs, F2Photon_IzNBars, f2photon_pts);
        double F2_chi2 = f2.chi2(F2_IzNs, F2_IzNBars, f2_pts);
        double FL_chi2 = fl.chi2(FL_IzNs, FL_IzNBars, fl_pts);
        double sigma_gp_chi2 = sigma_gp.chi2(sigma_gp_IzNs, sigma_gp_IzNBars, sigma_gp_pts);
        double sigma_gg_chi2 = sigma_gg.chi2(sigma_gg_IzNs, sigma_gg_IzNBars, sigma_gg_pts);
        double sigma_pp_chi2 = sigma_pp.chi2(sigma_pp_IzNs, sigma_pp_IzNBars, sigma_pp_pts);

        double chi2 = F2Photon_chi2 + F2_chi2 + FL_chi2 + sigma_gp_chi2 + sigma_gg_chi2 + sigma_pp_chi2;
        cout << "chi2: " << chi2 << endl; 
        return chi2;
    };


    // Start the fit now
    vector<double> X_guess = {k1, k2, k3, k4, kbar1, kbar2, kbar3, kbar4};
    double delta = 0.5;
    vector<double> X_opt = optimFunction(X_guess, f, delta);
    
    // Print X_opt
    cout << "Found a minimum for: ";
    for(int i = 0; i < X_opt.size(); i++) cout << X_opt[i] << " ";
    cout << endl;
    cout << "Number of experimental points: " << npoints << endl;
    cout << "Number of degrees of freedom: " << npoints - X_opt.size() << endl;
    cout << "Best chi2 / Ndof: " << f(X_opt) / (npoints - X_opt.size()) << endl;

    return 0;
}