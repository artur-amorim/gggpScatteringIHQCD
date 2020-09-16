#include <iostream>
#include <exception>
#include <string>
#include "gluonPDF.h"
#include "schrodinger/schrodinger.h"
#include "uncertainties.hpp"

using namespace std;

int main(int argc, char ** argv)
{  
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " nQ2s Q2 values hessian_path" << endl;
        throw runtime_error("Unsufficient number of parameters.");
    }
    // Create a structure with all the necessary kinematical points
    chebSetN(1000);
    const int nQ2s = stoi(argv[1]);
    const double log_xmin = -6;
    const double log_xmax = 0;
    const double h = 0.0001 ;

    vector<vector<double> > kinPts;
    vector<double> Q2s, xs;
    double Q2;
    for(int i = 1; i <= nQ2s; i ++)
    {
        Q2 = stod(argv[i+1]);
        for(double log_x = log_xmin; log_x <= log_xmax; log_x += h)
        {
            Q2s.push_back(Q2);
            xs.push_back(pow(10, log_x));
        }
    }
    kinPts = {Q2s, xs};

    GluonPDF gluonpdf;

    // Central value point: invls a b c d g1 g2 g3 g4
    vector<double> central_value = {6.46892, -4.69919, 1.12825, 0.664399, -0.0982592, -0.0509552, 0.0176742, -0.0739424, 0.35877};
    vector<double> xgs = gluonpdf.computePDF(kinPts, central_value);

    // Now we compute the uncertainties
     // Load the Hessian matrix
     arma::mat hessian;
     hessian.load(argv[2+nQ2s]);
     // Compute Mik
     arma::mat Mik = computeMik(hessian);
     // Define sqrtDeltaChi2 to 1
     const double sqrtDeltaChi2 = 10;

     // Now we compute the g(S^{+-}_k) and store them in a vector container
     vector<vector<double> > xgsVals = {xgs};
     // k is given by the number of columns of Mik
     const int n_pars = Mik.n_cols;

     vector<double> pars(central_value.size(), 0);
     // Iterate over k
     for(int k = 0; k < n_pars; k++)
     {
         // Start with S^{+}_k
         // compute a(S^{+}_k)_i = a0_i + \sqrt(DeltaChi2) M_{ik}
         for(int i = 0; i < central_value.size(); i++) pars[i] = central_value[i] + sqrtDeltaChi2 * Mik(i,k);
         // Compute g(S^{+}_k)
         xgsVals.push_back(gluonpdf.computePDF(kinPts, pars));
         // compute a(S^{-}_k)_i = a0_i - \sqrt(DeltaChi2) M_{ik}
         for(int i = 0; i < central_value.size(); i++) pars[i] = central_value[i] - sqrtDeltaChi2 * Mik(i,k);
         // Compute g(S^{-}_k)
         xgsVals.push_back(gluonpdf.computePDF(kinPts, pars));
     }

     //xgsVals = {central value, xg(S^{+}_1), xg(S^{-}_1), ...}

     // Compute \Delta xg  = 0.5 * sqrt(sum((xg(S+_k) - xg(S-_k))^2))
     vector<double> xgerror(kinPts[0].size(),0.0);
     double error;
     for(int i = 0; i < xgerror.size(); i++)
     {
         error = 0.0;
         for(int k = 1; k < xgsVals.size(); k += 2) error += pow(xgsVals[k][i] - xgsVals[k+1][i],2);
         error = 0.5 * sqrt(error);
         xgerror[i] = error;
     }

    ofstream myfile;
    string ihqcd_pdf_path;
    cout << "Insert file name where the IHQCD PDF predictions are saved." << endl;
    cin >> ihqcd_pdf_path;
    myfile.open(ihqcd_pdf_path);
    myfile << "Q2\tx\txg\txgError" << endl;
    for(int i = 0; i < kinPts[0].size(); i++)
        myfile << kinPts[0][i] << '\t' << kinPts[1][i] << '\t' << xgsVals[0][i] << '\t' << xgerror[i] << endl;
    myfile.close();
    
    return 0;
}