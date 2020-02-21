#include <iostream>
#include <fstream>
#include "IHQCD.h"
#include "U1NNMode.h"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/interpolation/Spline_Interp.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    double Q2 = stod(argv[1]);
    U1NNMode mode(Q2);
    cout << "Computing mode with Q2 = " << mode.Q2() << endl;
    mode.computeMode();
    cout << "Mode computed" << endl;
    mode.saveMode("Q2_"+to_string(mode.Q2())+".txt");

    return 0; 
}