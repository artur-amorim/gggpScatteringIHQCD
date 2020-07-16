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
#include "HardPomeron.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    // Kernel parameters
    double invls, a, b, c, d;
    // Gravitational photon couplings
    double k1, k2, k3, k4;
    // Gravitational proton couplings
    double kbar1, kbar2, kbar3, kbar4;
    string f2photon_path, f2_path, fl_path, sigma_gp_path, sigma_gg_path;
    if(argc < 19)
    {
        f2photon_path = "expdata/F2_photon/F2Photon_data.txt";
        f2_path = "expdata/DIS/F2_data.txt";
        fl_path = "expdata/DIS/FL_data.txt";
        sigma_gp_path = "expdata/GammaP/SigmaGammaP_PDG_data_W_gt_10.txt";
        sigma_gg_path = "expdata/GammaGamma/SigmaGammaGamma_PDG_data_W_gt_10.txt";
        invls = 6.491; a = -4.567; b = 1.485; c = 0.653; d = -0.113;
        k1 = 0.0482055; k2 = 0.0698707; k3 = 0.173638; k4 = -0.262936;
        kbar1 = -16.8249; kbar2 = -5.82798; kbar3 = 0.194873; kbar4 = -51.7903;
        cout << "Program usage: " + string(argv[0]) + " f2photon_path f2_path fl_path sigma_gp_path sigma_gg_path invls a b c d k1 k2 k3 k4 kbar1 kbar2 kbar3 kbar4" << endl;
        cout << "Using default values." << endl;
    }
    else
    {
        f2photon_path = argv[1]; f2_path = argv[2];
        fl_path = argv[3]; sigma_gp_path = argv[4];
        sigma_gg_path = argv[5];
        invls = stod(argv[6]); a = stod(argv[7]); b = stod(argv[8]); c = stod(argv[9]); d = stod(argv[10]);
        k1 = stod(argv[11]); k2 = stod(argv[12]); k3 = stod(argv[13]); k4 = stod(argv[14]);
        kbar1 = stod(argv[15]); kbar2 = stod(argv[16]); kbar3 = stod(argv[17]); kbar4 = stod(argv[18]);
    }
    cout << "Starting fit with values:" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << endl;
    cout << "Using the following gravitational couplings values" << endl;
    cout << "k1: " << k1 << " k2: " << k2 << " k3: " << k3 << " k4: " << k4 << endl;
    cout << "kbar1: " << kbar1 << " kbar2: " << kbar2 << " kbar3: " << kbar3 << " kbar4: " << kbar4 << endl;

    // Define the process observables and load the data needed for the fit
    F2Photon f2photon(f2photon_path);
    F2 f2(f2_path); FL fl(fl_path);
    SigmaGammaProton sigma_gp(sigma_gp_path);
    SigmaGammaGamma sigma_gg(sigma_gg_path);

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

    // Setup Chebyschev computation
    chebSetN(1000);

    // Setup HardPomeron Kernel and compute the Reggeons
    HardPomeron hard(4, {invls, a, b, c, d});
    const vector<double> ks = {k1, k2, k3, k4};
    const vector<double> kbars = {kbar1, kbar2, kbar3, kbar4};


    auto f = [&f2photon, &f2, &fl, &sigma_gp, &sigma_gg, &f2photon_pts, &f2_pts, &fl_pts, &sigma_gp_pts, &sigma_gg_pts, &hard, &ks, &kbars] (const std::vector<double> & X)
    {
        // X is a vector with the values of invls a b c d k1, k2, k3, k4, kbar1, kbar2, kbar3, kbar4
        vector<double> kernelPars = {X[0], X[1], X[2], X[3], X[4]};

        cout << "invls: " << X[0] << " a: " << X[1] << " b: " << X[2] << " c: " << X[3] << " d: " << X[4] << endl;

        hard.computeReggeTrajectories(kernelPars);
        vector<Reggeon> reggeons = computeReggeons(hard, 0, 4);
        Spectra spec(0, reggeons);
        vector<Spectra> spectrum = {spec};

        // Compute IzNs
        vector<kinStruct> F2Photon_IzNs = f2photon.getIzs(f2photon_pts, spectrum);
        vector<kinStruct> F2_IzNs = f2.getIzs(f2_pts, spectrum);
        vector<kinStruct> FL_IzNs = fl.getIzs(fl_pts, spectrum);
        vector<kinStruct> sigma_gp_IzNs = sigma_gp.getIzs(sigma_gp_pts, spectrum);
        vector<kinStruct> sigma_gg_IzNs = sigma_gg.getIzs(sigma_gg_pts, spectrum);

        // Now we need to compute the gn's according to the definition in the notes
        
        vector<double> Im_gn_gg(reggeons.size()), Im_gn_gp(reggeons.size());
        for(int i = 0; i < reggeons.size(); i++)
        {
            double jn = reggeons[i].getJ(), djndt = reggeons[i].getdJdt();
            // Common factor
            complex<double> gn(1/tan(M_PI_2 * jn), 1);
            gn = - M_PI_2 * gn * djndt / pow(2, jn);
            // Now we make the specific computations
            Im_gn_gg[i] = imag(gn * ks[i] * ks[i] * sigma_gg_IzNs[0].izns[i]);
            Im_gn_gp[i] = imag(gn * ks[i] * kbars[i]);
            if(i == 0 or i == 3)
            {
                Im_gn_gg[i] = - Im_gn_gg[i]; Im_gn_gp[i] = - Im_gn_gp[i];
            }
        }

        // Compute the IzNsBar
        vector<kinStruct> F2Photon_IzNBars = f2photon.getIzsBar(f2photon_pts, spectrum, Im_gn_gg);
        vector<kinStruct> F2_IzNBars = f2.getIzsBar(f2_pts, spectrum, Im_gn_gp);
        vector<kinStruct> FL_IzNBars = fl.getIzsBar(fl_pts, spectrum, Im_gn_gp);
        vector<kinStruct> sigma_gp_IzNBars = sigma_gp.getIzsBar(sigma_gp_pts, spectrum, Im_gn_gp);
        vector<kinStruct> sigma_gg_IzNBars = sigma_gg.getIzsBar(sigma_gg_pts, spectrum, Im_gn_gg);

        // Compute the chi2 of each process
        double F2Photon_chi2 = f2photon.chi2(F2Photon_IzNs, F2Photon_IzNBars, f2photon_pts);
        double F2_chi2 = f2.chi2(F2_IzNs, F2_IzNBars, f2_pts);
        double FL_chi2 = fl.chi2(FL_IzNs, FL_IzNBars, fl_pts);
        double sigma_gp_chi2 = sigma_gp.chi2(sigma_gp_IzNs, sigma_gp_IzNBars, sigma_gp_pts);
        double sigma_gg_chi2 = sigma_gg.chi2(sigma_gg_IzNs, sigma_gg_IzNBars, sigma_gg_pts);

        double chi2 = F2Photon_chi2 + F2_chi2 + FL_chi2 + sigma_gp_chi2 + sigma_gg_chi2;
        cout << "chi2: " << chi2 << endl; 
        return chi2;
    };

    // Start the fit now
    vector<double> X_guess = {invls, a, b, c, d};
    vector<double> delta = {0.5, 0.5, 0.5, 0.5, 0.5};
    vector<double> X_opt = optimFunction(X_guess, f, delta);
    
    // Print X_opt
    cout << "Found a minimum for: ";
    for(int i = 0; i < X_opt.size(); i++) cout << X_opt[i] << " ";
    cout << endl;
    cout << "Number of experimental points: " << npoints << endl;
    cout << "Number of degrees of freedom: " << npoints - X_opt.size() << endl;
    cout << "Best chi2 / Ndof: " << f(X_opt) / (npoints - X_opt.size()) << endl;

    hard.computeReggeTrajectories({X_opt[0], X_opt[1], X_opt[2], X_opt[3], X_opt[4]});
    vector<Reggeon> reggeons = computeReggeons(hard, 0, 4);
    Spectra spec(0, reggeons);
    vector<Spectra> spectrum = {spec};
    
    vector<kinStruct> sigma_gg_IzNs = sigma_gg.getIzs(sigma_gg_pts, spectrum);

    cout << "Values of the imaginary parts:" << endl;

    vector<double> Im_gn_gg(reggeons.size()), Im_gn_gp(reggeons.size());
    for(int i = 0; i < reggeons.size(); i++)
    {
        double jn = reggeons[i].getJ(), djndt = reggeons[i].getdJdt();
        // Common factor
        complex<double> gn(1/tan(M_PI_2 * jn), 1);
        gn = - M_PI_2 * gn * djndt / pow(2, jn);
        // Now we make the specific computations
        Im_gn_gg[i] = imag(gn * ks[i] * ks[i] * sigma_gg_IzNs[0].izns[i]);
        Im_gn_gp[i] = imag(gn * ks[i] * kbars[i]);
        if(i == 0 or i == 3)
        {
            Im_gn_gg[i] = - Im_gn_gg[i]; Im_gn_gp[i] = - Im_gn_gp[i];
        }
    }
    cout << "Im_gn_gg:" << endl;
    for(int i = 0; i < Im_gn_gg.size(); i++) cout << Im_gn_gg[i] << '\t';
    cout << endl;
    cout << "Im_gn_gp:" << endl;
    for(int i = 0; i < Im_gn_gp.size(); i++) cout << Im_gn_gp[i] << '\t';
    cout << endl;

    return 0;
}