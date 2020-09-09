#include <iostream>
#include <string>
#include <complex>
#include <random>
#include "IHQCD.h"
#include "F2.h"
#include "FL.h"
#include "SigmaGammaProton.h"
#include "HardPomeron.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"
#include "uncertainties.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    arma::mat hessian;
    if(argc < 2)
    {
        cout << "Program usage: " + string(argv[0]) + " f2_path fl_path sigma_gp_path invls a b c d Im_g0 Im_g1 Im_g2 Im_g3 h" << endl;
        cout << "Or " << string(argv[0]) << " hessian_path" << endl;
        return 0;
    }
    else if (argc == 2)
    {
        hessian.load(string(argv[1]));
    }
    else if (argc == 14)
    {
        string f2_path = argv[1], fl_path = argv[2], sigma_gp_path = argv[3];
        // Gluon Kernel parameters
        double invls = stod(argv[4]), a = stod(argv[5]), b = stod(argv[6]), c = stod(argv[7]), d = stod(argv[8]);
        // Gravitational photon couplings
        double Im_g0 = stod(argv[9]), Im_g1 = stod(argv[10]), Im_g2 = stod(argv[11]), Im_g3 = stod(argv[12]);
        double h = stod(argv[13]);
    
        cout << "Computing uncertainties for values:" << endl;
        cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
        cout << "Im_g0: " << Im_g0 << " Im_g1: " << Im_g1 << " Im_g2: " << Im_g2 << " Im_g3: " << Im_g3 << endl;

        // Define the process observables and load the data needed for the fit
        F2 f2(f2_path);
        FL fl(fl_path);
        SigmaGammaProton sigma_gp(sigma_gp_path);

        // Get the experimental points
        vector<vector<double> > f2_pts = f2.expKinematics(), fl_pts = fl.expKinematics();
        vector<vector<double> > sigma_gp_pts = sigma_gp.expKinematics();
    
        // Setup Chebyschev computation
        chebSetN(1000);

        // Setup HardPomeron Kernel and compute the Reggeons
        HardPomeron hard(4 , {invls, a, b, c, d});


        // Definition of the functions we want to fit
        auto f = [&f2, &fl, &sigma_gp, &f2_pts, &fl_pts, &sigma_gp_pts, &hard] (const std::vector<double> & X)
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
            vector<kinStruct> sigma_gp_IzNs = sigma_gp.getIzs(sigma_gp_pts, spectrum);

            // Compute the IzNsBar
            vector<kinStruct> F2_IzNBars = f2.getIzsBar(f2_pts, spectrum, Im_gns);
            vector<kinStruct> FL_IzNBars = fl.getIzsBar(fl_pts, spectrum, Im_gns);
            vector<kinStruct> sigma_gp_IzNBars = sigma_gp.getIzsBar(sigma_gp_pts, spectrum, Im_gns);

            // Compute the chi2 of each process
            double F2_chi2 = f2.chi2(F2_IzNs, F2_IzNBars, f2_pts);
            double FL_chi2 = fl.chi2(FL_IzNs, FL_IzNBars, fl_pts);
            double sigma_gp_chi2 = sigma_gp.chi2(sigma_gp_IzNs, sigma_gp_IzNBars, sigma_gp_pts);
            double chi2 = F2_chi2 + FL_chi2 + sigma_gp_chi2;
            cout << "chi2: " << chi2 << endl; 
            return chi2;
        };

        // Compute the hessian matrix and save it
        vector<double> X = {invls, a, b, c, d, Im_g0, Im_g1, Im_g2, Im_g3};
        hessian = computeHessian(f, X, h);
        string hessian_path;
        cout << "Introduce the file name that will store the hessian of gamma proton processes." << endl;
        cin >> hessian_path;
        saveMatrix(hessian, hessian_path);
    }
    else
    {
        cout << "Invalid number of arguments. Call this program without arguments for instructions." << endl;
        return 0;
    }
    
    // Compute the parameter uncertainties
    vector<double> parUncert = parameterUncertainties(hessian);

    // Print the parameter values and the corresponding uncertainties
    cout << "Delta_invls: " << parUncert[0] << endl;
    cout << "Delta_a: " << parUncert[1] << endl;
    cout << "Delta_b: " << parUncert[2] << endl;
    cout << "Delta_c: " << parUncert[3] << endl;
    cout << "Delta_d: " << parUncert[4] << endl;
    cout << "Delta_Im_g0: " << parUncert[5] << endl;
    cout << "Delta_Im_g1: " << parUncert[6] << endl;
    cout << "Delta_Im_g2: " << parUncert[7] << endl;
    cout << "Delta_Im_g3: " << parUncert[8] << endl;

    return 0;
}