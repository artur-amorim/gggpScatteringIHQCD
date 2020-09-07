#include <string>
#include <vector>
#include <fstream>
#include <functional>
#include <cmath>
#include <boost/math/distributions/chi_squared.hpp>
#include <armadillo>

#include "schrodinger/common.h"

template <class T>
arma::mat computeHessian(T &func, const std::vector<double> &point, const double h)
{
    /*
        Computes numerically the hessian matrix of the function func at point.
        func is a multivariable function that takes as input a std::vector<double>
        point is a std::vector<double> where we want to evaluate the Hessian matrix of func
        The routine starts to compute the dimension of the Hessian matrix. This is done
        by knowing the size of point. The Hessian components (partial derivatives of func)
        are computed using the central difference approximation (Abramowitz and Stegun 1972, p. 884)
        d2f/dxi^2 = (-f(x+2 hi ei) + 16f(x + hi ei) - 30 f(x) + 16 f(x - hi ei) - f(x - 2 hi ei)) / (12 hi^2)
        d2f/(dxjdxi) = (f(x + hi ei + hj ej) - f(x + hi ei - hj ej) - f(x - hi ei + hj ej) + f(x - hi ei - hj ej))/(4 * hi hj)
    */
   const int n = point.size();
   arma::mat hessian(n,n, arma::fill::zeros);
   // Variable that will store the value of func at point. Usefull if evaluating func is expensive
   double fpt = func(point);
   // Auxiliary points to compute Hessian elements
   std::vector<double> p1, p2, p3, p4;
   for(int i = 0; i < n; i++)
   {
       for(int j = 0; j < n; j++)
       {
           if(i == j)
           {
               // Diagonal element
               p1 = point; p2 = point; p3 = point; p4 = point;
               p1[i] += 2 * h; p2[i] += h; p3[i] -= h; p4[i] -= 2*h;
               hessian(i,i) = (-func(p1) + 16.0 * func(p2) - 30. * fpt + 16.0 * func(p3) - func(p4)) / (12 * h * h) ;
           }
           else if ( j > i)
           {
               // upper off-diagonal element
               p1 = point; p2 = point; p3 = point; p4 = point;
               p1[i] += h; p1[j] += h; p2[i] += h; p2[j] -= h; p3[i] -= h; p3[j] += h; p4[i] -= h; p4[j] -= h;
               hessian(i, j) = 0.25 * (func(p1) - func(p2) - func(p3) + func(p4)) / (h*h) ;
           }
           else hessian(i, j) = hessian(j, i);                  // The Hessian matrix is symmetric
       }
   }
   return hessian;
}

arma::mat computeCovarianceMatrix(const arma::mat &hessian)
{
    /*
        Given the Hessian matrix H, it computes the covariant matric C = \alpha^{-1},
        where \alpha = 0.5 * H
    */
   arma::mat alpha = 0.5 * hessian;
   arma::mat C = arma::inv(alpha);
   return C;
}

void saveMatrix(const arma::mat &matrix, const std::string &filepath)
{
    /*
        Saves the matrix in a file stored in filepath
    */
   const int n_rows = matrix.n_rows;
   const int n_cols = matrix.n_cols;
   std::ofstream myfile;
   myfile.open(filepath);
   for(int i = 0; i < n_rows; i++)
   {
       for(int j = 0; j < n_cols; j++) myfile << matrix(i, j) << '\t';
       myfile << std::endl;
   }
   myfile.close();
}

arma::mat computeMik(const arma::mat &hessian)
{
    /*
        Given the Hessian matrix it computes the matrix of rescaled eigenvectors Mik.
        The routine has the following steps:
            - Computation of the eigenvalues and eigenvectors of hessian through arma::eig_sym (Hessian is a symmetric matrix)
            - The eigenvectors will be columns of the matrix eigvec defined below
            - Rescale the eigenvectors by \sqrt(2/\lambda) where \lambda is the corresponding eigenvector
            - return the rescaled matrix.
    */
   const int n_pars = hessian.n_rows;
   // Compute the eigensystem of the Hessian matrix
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, hessian);

    // Rescale the eigenvectors. eigvec will become M_{ik}
    for(int i = 0; i < n_pars; i++) eigvec.col(i) *= std::sqrt(2.0/eigval(i)) ;

    return eigvec;
}

double computeDeltaChi2(int n_deg_freed, const double conf_limit)
{
    // Computes DeltaChi2 given the number of parameters and confidence limit
    boost::math::chi_squared_distribution<double> dist(n_deg_freed);
    std::function<double(double)> f = [dist, conf_limit] (const double &x)
    {
        return cdf(dist, x) - conf_limit;
    };
    double Deltachi2 = zbrent(f, 0, 100, 1e-9);
    return Deltachi2;
}

std::vector<double> parameterUncertainties(const arma::mat &hessian)
{
    /*
        Given the Hessian matrix it compute the parameter uncertainties of a fit.
        The parameter uncertainties are the ones defined by the boundary 
        of the d-dimensional region in the parameter space such that chi2_fit + 1.
        It assumes that around the minimum the chi2 function is well approximated quadratically.
        The routine starts by defining \Delta chi2 = 1.
        Then we compute the eigenvalues and orthonormal eigenvectors of the hessian.
        The kth eigenvector is then rescaled by sk = \sqrt( 2 / \eps_k ) in such
        a way that \Delta chi2 <= 0.5 H_{ij}(a_i - a0_i)(a_j - a0_j) = z_i z_i 
        where an implicit sum is assmued.
        The z_i are the coordinates of a point in the parameter space relative to
        the basis of eigenvectors of the hessian matrix.
        A matrix M_{ik} is constructed from the rescaled eigenvectors having them as columns. M_{ik}
        is the ith component of the kth eigenvector in the canonical basis.
        The parameter uncertainties are then given by \Delta a_i = \sqrt(\Delta chi2)\sqrt(Mik Mik)
        where a summation over k is implicit.
        A vector containing the \Delta a_i is returned.
    */
    // Computation of DeltaChi2
    const int n_pars = hessian.n_rows;
    double Deltachi2 = 1;

    // Computation of the matrix Mik
    arma::mat Mik = computeMik(hessian);

    // \Delta a_i = \sqrt(\Delta chi2)\sqrt(Mik Mik)
    std::vector<double> parUncert(n_pars,0.0);
    for(int i = 0; i < n_pars; i++)
    {
        for(int j = 0; j < n_pars; j++) parUncert[i] += Mik(i,j) * Mik(i, j);
        parUncert[i] = std::sqrt(Deltachi2) * std::sqrt(parUncert[i]);
    }
    return parUncert;
}