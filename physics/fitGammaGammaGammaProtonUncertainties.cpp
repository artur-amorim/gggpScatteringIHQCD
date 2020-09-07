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
#include "uncertainties.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    arma::mat hessian;
    if(argc < 2)
    {
        cout << "Program usage: " + string(argv[0]) + " f2photon_path f2_path fl_path sigma_gp_path sigma_gg_path k1 k2 k3 k4 kbar1 kbar2 kbar3 kbar4 h" << endl;
        cout << "Or " << argv[0] << " hessian_path" << endl;
        return 0;
    }
    else if (argc == 2)
    {
        hessian.load(string(argv[1])); 
    }
    else if (argc == 15)
    {
        string f2photon_path = argv[1], f2_path = argv[2], fl_path = argv[3], sigma_gp_path = argv[4], sigma_gg_path = argv[5];
        // Gravitational photon couplings
        double k1 = stod(argv[6]), k2 = stod(argv[7]), k3 = stod(argv[8]), k4 = stod(argv[9]);
        // Gravitational proton couplings
        double kbar1 = stod(argv[10]), kbar2 = stod(argv[11]), kbar3 = stod(argv[12]), kbar4 = stod(argv[13]);
        double h = stod(argv[14]);
        
        cout << "Computing uncertainties for values:" << endl;
        cout << "k1: " << k1 << " k2: " << k2 << " k3: " << k3 << " k4: " << k4 << endl;
        cout << "kbar1: " << kbar1 << " kbar2: " << kbar2 << " kbar3: " << kbar3 << " kbar4: " << kbar4 << endl;

        // Define the process observables and load the data needed for the fit
        F2Photon f2photon(f2photon_path);
        F2 f2(f2_path); FL fl(fl_path);
        SigmaGammaProton sigma_gp(sigma_gp_path);
        SigmaGammaGamma sigma_gg(sigma_gg_path);

        // Get the experimental points
        vector<vector<double> > f2photon_pts = f2photon.expKinematics();
        vector<vector<double> > f2_pts = f2.expKinematics() , fl_pts = fl.expKinematics();
        vector<vector<double> > sigma_gp_pts = sigma_gp.expKinematics(), sigma_gg_pts = sigma_gg.expKinematics();

        // Setup Chebyschev computation
        chebSetN(1000);

        // Setup HardPomeron Kernel and compute the Reggeons
        HardPomeron hard(4, {6.46892, -4.69919, 1.12825, 0.664399, -0.0982592});
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

        auto f = [&f2photon, &f2, &fl, &sigma_gp, &sigma_gg, &spectrum, &f2photon_pts, &f2_pts, &fl_pts, &sigma_gp_pts, &sigma_gg_pts,
             &F2Photon_IzNs, &F2_IzNs, &FL_IzNs, &sigma_gp_IzNs, &sigma_gg_IzNs, &reggeons] (const std::vector<double> & X)
        {
            // X is a vector with the values of k1, k2, k3, k4, kbar1, kbar2, kbar3, kbar4
            vector<double> ks = {X[0], X[1], X[2], X[3]};
            vector<double> kbars = {X[4], X[5], X[6], X[7]};
            /*
                Now we need to compute the gn's according to the definition in the notes
            */
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
            return chi2;
        };

        // Compute the uncertainties
        vector<double> X = {k1, k2, k3, k4, kbar1, kbar2, kbar3, kbar4};
        hessian = computeHessian(f, X, h);
        saveMatrix(hessian, "gamma_gamma_gamma_proton_hessian.txt");
    }
    else
    {
        cout << "Invalid number of arguments. Call this program without arguments for instructions." << endl;
        return 0;
    }

    // We can now compute the parameter uncertainties
    vector<double> parUncert = parameterUncertainties(hessian);

    // Print the parameter values and the corresponding uncertainties
    cout << "Delta_k1: " << parUncert[0] << endl;
    cout << "Delta_k2: " << parUncert[1] << endl;
    cout << "Delta_k3: " << parUncert[2] << endl;
    cout << "Delta_k4: " << parUncert[3] << endl;
    cout << "Delta_kbar1: " << parUncert[4] << endl;
    cout << "Delta_kbar2: " << parUncert[5] << endl;
    cout << "Delta_kbar3: " << parUncert[6] << endl;
    cout << "Delta_kbar4: " << parUncert[7] << endl;

    return 0;
}