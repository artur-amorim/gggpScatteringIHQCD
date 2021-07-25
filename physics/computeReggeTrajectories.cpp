#include <iostream>
#include <fstream>
#include <cmath>
#include "HardPomeron.h"
#include "IHQCD.h"
#include "schrodinger/schrodinger.h"
#include "methods/interpolation/Spline_Interp.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    double invls, a, b, c, d;
    int N;
    string regge_traj_path;
    if (argc < 8)
    {
        invls = 6.47; a = -4.70; b = 1.13; c =  0.66; d = -0.098;
        N = 1000;
        regge_traj_path = "regge_trajectories";
    }
    else
    {
        invls = stod(argv[1]); a = stod(argv[2]); b = stod(argv[3]); c = stod(argv[4]); d = stod(argv[5]);
        N = stoi(argv[6]); regge_traj_path = argv[7];
    }

    cout << "Computing the gluon Regge trajectories with:" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
    
    vector<double> gluon_pars = {invls, a, b, c, d};
    HardPomeron gluon(4, gluon_pars);

    // initialise Chebyschev matrices
    chebSetN(N);

    // Compute the trajectories 
    cout << "Computing the trajectories." << endl;
    vector<double> z = ihqcd().getZ(), Vgluon(0);
    vector<double> ts(0);
    
    ofstream gluon_traj;
    gluon_traj.open(regge_traj_path + "_gluon.txt");
    gluon_traj << "j\tt1(j)\tt2(j)\tt3(j)\tt4(j)" << endl;

    vector<double> js(0);
    vector<vector<double>> ts_gluon(4, vector<double>(0));

    for(double j = 0.1; j <= 4.1; j += 0.01)
    {
        // Append value of j
        js.push_back(j);
        // Compute the potentials
        Vgluon = gluon.computePotential(j);

        // Compute the spectrum of the gluon
        List spec = computeSpectrum(z, Vgluon, 4);
        ts = spec.Es;
        gluon_traj << j << '\t' << ts[0] << '\t' << ts[1] << '\t' << ts[2] << '\t' << ts[3] << endl;
        // Append the t values to ts_gluon
        ts_gluon[0].push_back(ts[0]); ts_gluon[1].push_back(ts[1]); 
        ts_gluon[2].push_back(ts[2]); ts_gluon[3].push_back(ts[3]);
        
    }
    gluon_traj.close();

    // Compute the intercepts of the Pomeron and meson trajectories
    Spline_Interp<double> gluon_traj_0(ts_gluon[0], js), gluon_traj_1(ts_gluon[1], js), gluon_traj_2(ts_gluon[2], js), gluon_traj_3(ts_gluon[3], js) ;

    const double j0g = gluon_traj_0.interp(0), j1g = gluon_traj_1.interp(0), j2g = gluon_traj_2.interp(0), j3g = gluon_traj_3.interp(0);

    cout << "Gluon Intercepts:" << endl;
    cout << "j0g = " << j0g << " j1g = " << j1g << " j2g = " << j2g << " j3g = " << j3g << endl;

    return 0;
}