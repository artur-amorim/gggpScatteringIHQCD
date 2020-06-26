#include <iostream>
#include <cmath>
#include "IHQCD.h"
#include "methods/vectorOperators.hpp"
#include "schrodinger/schrodinger.h"

using namespace std;

int main()
{
    vector<double> z = ihqcd().getZ();
    vector<double> e2As = ihqcd().getE2As();
    vector<double> dAstring = ihqcd().getdAstring(), d2Astring = ihqcd().getd2Astring();
    vector<double> dPhi = ihqcd().getdPhi(), d2Phi = ihqcd().getd2Phi();
    vector<double> VSG = ihqcd().getU0();
    vector<double> VTG = ihqcd().getU2();
    vector<double> Vf2 = 0.5*(3.0*d2Astring-d2Phi) + 0.25*(3.0*dAstring-dPhi)*(3.0*dAstring-dPhi) + 4.0 * e2As;

    chebSetN(1000);
    List spectrum_SG = computeSpectrum(z, VSG, 4, "cheb");
    List spectrum_TG = computeSpectrum(z, VTG, 4, "cheb");
    List spectrum_f2 = computeSpectrum(z, Vf2, 4, "cheb");
    vector<double> m2_SG = spectrum_SG.Es, m2_TG = spectrum_TG.Es, m2_f2 = spectrum_f2.Es;
    cout << "SCALAR GLUEBALL SPECTRUM:" << endl;
    for(int i = 0; i < m2_SG.size(); i++) cout << sqrt(m2_SG[i]) << '\t';
    cout << endl;
    cout << "TENSOR GLUEBALL SPECTRUM:" << endl;
    for(int i = 0; i < m2_TG.size(); i++) cout << sqrt(m2_TG[i]) << '\t';
    cout << endl;
    cout << "TENSOR MESON SPECTRUM:" << endl;
    for(int i = 0; i < m2_f2.size(); i++) cout << sqrt(m2_f2[i]) << '\t';
    cout << endl;

    return 0;
}
