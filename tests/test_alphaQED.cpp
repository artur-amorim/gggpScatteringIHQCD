#include <iostream>
#include <fstream>
#include <cmath>
#include "alphaQED.h"

using namespace std;

int main()
{
    // We are interested in the mZ^2 scale
    double s = 91.1876 * 91.1876;
    double alpha = alphaQED(s);
    double Invalpha = InvalphaQED(s);
    cout << "alphaQED(mZ^2): " << alpha << endl;
    cout << "InvAlphaQED(mZ^2): " << Invalpha << endl;

    cout << "InvAlphaQED(-2.1 GeV^2) - InvAlphaQED(-6.2 GeV^2: " << InvalphaQED(2.1) - InvalphaQED(6.2) << endl;
    cout << "InvAlphaQED(-12.25 GeV^2) - InvAlphaQED(-3434 GeV^2): " << InvalphaQED(12.25) - InvalphaQED(3434) << endl;  

    cout << endl;
    cout << "Computing alphaQED for other values of Q2. They will be saved in the file alphaQED.txt" << endl;

    ofstream myfile;
    myfile.open("alphaQED.txt");
    
    myfile << "Q2\talphaQED\tInvalphaQED" << endl;
    for(double log10Q2 = 0; log10Q2 <= 4; log10Q2 += 0.01)
    {
        double Q2 = pow(10, log10Q2);
        myfile << Q2 << '\t' << alphaQED(Q2) << '\t' << InvalphaQED(Q2) << endl;
    }
    myfile.close();

    return 0;
}