#include <iostream>
#include <fstream>
#include <vector>
#include "Reggeon.h"
#include "GluonKernel.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main()
{
    // Compute Chebyschev matrices
    chebSetN(800);

    // Setup gluon kernel and compute the Reggeons for t = 0
    vector<double> gluon_pars = {1/0.178};
    GluonKernel gluon(2, gluon_pars);
    gluon.computeReggeTrajectories();
    vector<Reggeon> reggeons = computeReggeons(gluon, 0.0, 2);
    // Print the J of each reggeon and save the wave function
    for(int i = 0; i < reggeons.size(); i++)
    {
        Reggeon reg = reggeons[i];
        cout << reg.getName() << '\t' << reg.getJ() << '\t' << reg.getdJdt() << endl;
        vector<vector<double> > wf = reg.getWf();
        vector<double> z = wf[0], psi = wf[1];
        ofstream file;
        file.open("wf_"+to_string(reg.getIndex())+".txt");
        file << "z" << '\t' << "psi" << endl;
        for(int j = 0; j < z.size(); j++) file << z[j] << '\t' << psi[j] << endl;
        file.close();
    }

    return 0;
}