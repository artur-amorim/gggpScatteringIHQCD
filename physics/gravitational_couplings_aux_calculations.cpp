#include <iostream>
#include <vector>
#include "Reggeon.h"
#include "HardPomeron.h"
#include "SigmaGammaGamma.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double invls, a, b, c, d;
    if (argc < 6)
    {
        invls = 6.491; a = -4.567; b = 1.485; c = 0.653; d = -0.113;
    }
    else
    {
        invls = stod(argv[1]);
        a = stod(argv[2]); b = stod(argv[3]); c = stod(argv[4]); d = stod(argv[5]);
    }
    chebSetN(1000);
    cout << "Computing the HardPomeorn kernel with following parameters" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;  

    // Setup HardPomeron Kernel
    HardPomeron hard(4, {invls, a, b, c, d});
    hard.computeReggeTrajectories();

    vector<Reggeon> reggeons = computeReggeons(hard, 0, 4);

    cout << "n\tjn\tdjn/dt" << endl;
    for(int i = 0; i < reggeons.size(); i++) cout << i << '\t' << reggeons[i].getJ() << '\t' << reggeons[i].getdJdt() << endl;

    // We now compute the IzN integrals
    SigmaGammaGamma sigma_gg("expdata/GammaGamma/gammagammaScatteringL3Processed.txt");
    cout << "n\tIzN" << endl;
    for(int i = 0; i < reggeons.size(); i++) cout << i << '\t' << sigma_gg.IzN({}, reggeons[i]) << endl;


    return 0;
}