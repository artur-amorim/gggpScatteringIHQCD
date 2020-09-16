#include <iostream>
#include <fstream>
#include <vector>
#include "Reggeon.h"
#include "HardPomeron.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main()
{
    // Compute Chebyschev matrices
    chebSetN(1000);

    // Setup gluon kernel and compute the Reggeons for t = 0
    HardPomeron hard(4, {6.46892, -4.69919, 1.12825, 0.664399, -0.0982592});
    hard.computeReggeTrajectories();
    vector<Reggeon> reggeons = computeReggeons(hard, 0.0, 4);
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