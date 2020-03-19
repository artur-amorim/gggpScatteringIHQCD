#include <iostream>
#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>

#include "IHQCD.h"
#include "alphaQED.h"
#include "F2Photon.h"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/vectorOperators.hpp"
#include "methods/search.hpp"

extern"C"
{
    void dqags_(double (*f)(double *, void *), void * params,
                double * a, double * b, double * epsabs, double * epsrel, 
                double * result, double * abserr, int * neval,int * ier,
                int * limit, int * lenw, int * last, int * iwork, double * work);
}

class F2PhotonIzNIntegrand
{
    private:
        Poly_Interp<double> f1, f3;
        U1NNMode f2;
    public:
        F2PhotonIzNIntegrand(const Poly_Interp<double> &func1, const U1NNMode &func2, const Poly_Interp<double> &func3);
        double operator()(const double x);
};

F2PhotonIzNIntegrand::F2PhotonIzNIntegrand(const Poly_Interp<double> &func1, const U1NNMode &func2, const Poly_Interp<double> &func3): f1(func1), f2(func2), f3(func3) {}

double F2PhotonIzNIntegrand::operator()(const double x) {return f1.interp(x) * f2.factor(x) * f3.interp(x);}

double fF2Photon(double * x, void * params)
{
    return ((F2PhotonIzNIntegrand *) params)->operator()(*x);
}

void F2Photon::computeU1NNModes()
{
    /*
        Computes the relevant U(1) Nonnormalizable modes from the know values of Q2
        present in the data points.
        The modes are stored in the std::vector container modes and are ordered by Q2.
    */
    // First clear the modes container
    modes.clear();
    // Get a std::vector with a unique list of Q2
    std::vector<std::vector<double> > data = this->getDataPts();
    std::vector<double> Q2s(0);
    if(data.size()!= 0) Q2s = data[0];
    else return ;
    std::sort(Q2s.begin(), Q2s.end());
    Q2s.erase(unique(Q2s.begin(), Q2s.end()), Q2s.end());
    // For each Q2 compute a mode
    std::cout << "Computing the U(1) Nonnormalizable modes of DIS object:" << std::endl;
    double Q2;
    for(int i = 0; i < Q2s.size(); i++)
    {
        Q2 = Q2s[i];
        std::cout << "Computing mode for Q2 = " << Q2 << std::endl;
        // Compute the mode
        U1NNMode mode(Q2);
        mode.computeMode();
        modes.push_back(mode);
    }
}

U1NNMode F2Photon::searchMode(const double Q2)
{
    /*
        Function that searches for the mode with virtuality Q2.
        Uses the binary search algorithm to find the mode.
        If the mode is not found it is computed.
    */
   U1NNMode mode(Q2);
   try
   {
       return binary_search<U1NNMode>(modes, mode);
   }
   catch(const char *)
   {
       std::cout << Q2 << std::endl;
       mode.computeMode();
       return mode;
   }
}

void F2Photon::loadData(std::string file_path)
{
     // Vector container with Q2, x, Fi and Fierr values, with i = 2, L
    if (file_path.size() == 0) return;
    std::ifstream file;
    file.open(file_path);
    if(file.fail()) std::runtime_error("Error: File with DIS data not opened."); 
    std::string line ;
    getline(file, line) ;
    std::cout << "Loading DIS data" << std::endl;
    std::vector<double> Q2, x, Fi, Fierr;
    std::vector<std::string> result;
    while(getline(file, line))
    {
        boost::split(result, line, boost::is_any_of("\t") ) ;
        x.push_back(stod(result[0])) ;
        Q2.push_back(stod(result[3])) ;
        Fi.push_back(stod(result[6])) ;
        Fierr.push_back(stod(result[7]));
    }
    // Check that all std::vector containers have the same side
    if( Fi.size() != Q2.size() || Q2.size() != x.size() || x.size() != Fierr.size())
    {
        std::cout << "Vector containers don't have the same size" << std::endl;
        throw "error";
    }
    this->setDataPts({Q2, x, Fi, Fierr}) ;
    std::cout << "Loaded F2Photon data." << std::endl;
}

void F2Photon::copy(const F2Photon &rhs)
{
    PhotonScattering::copy(rhs);
    modes = rhs.modes;
}

F2Photon::F2Photon(std::string file_path):
    PhotonScattering(), modes({})
{
    // Load data
    loadData(file_path);
    // Compute all the necessary U(1) NN modes
    // First setup its computation
    computeU1NNModes();
}

F2Photon::F2Photon(const F2Photon &f2):
    PhotonScattering(f2), modes(f2.modes){}

std::vector<double> F2Photon::expVal()
{
    // F2 data is present in index 2 of datapts
    return this->getDataPts()[2] ;
}

std::vector<double> F2Photon::expErr()
{
    // F2 err is present in index 3 of datapts
    return this->getDataPts()[3] ;
}

std::vector < std::vector<double> > F2Photon::expKinematics()
{
    std::vector < std::vector < double > > ans ;
    // Get the Q2s std::vector
    ans.push_back(this->getDataPts()[0]) ;
    // Get the xs std::vector
    ans.push_back(this->getDataPts()[1]) ;
    return ans ;
}

double F2Photon::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    // Get the kinematical values and J
    const double Q2 = kin[0];
    const double J = reg.getJ();
    std::vector<std::vector<double> > wf = reg.getWf();
    // Search for the correct mode
    std::vector<double> z = ihqcd().getZ();
    std::vector<double> fact1 = exp((1.5-J) * ihqcd().getAstring());
    Poly_Interp<double> f1(z, fact1, 4);
    U1NNMode mode = searchMode(Q2);
    // Interpolate the wavefunction
    Poly_Interp<double> f3(wf[0], wf[1], 4);
    // Define the integrand object
    F2PhotonIzNIntegrand integrand(f1, mode, f3);
    void * params = &integrand;
    // Compute the integral
    double izn = 0.0, abserr = 0.0;
    double a = z[0], b = z.back();
    // Absolute and relative tolerances desidered
    double epsabs = 1e-9, epsrel = 1e-9;
    if (Q2 < 2) epsabs = 1e-7;
    // Output variables 
    int neval = 0, ier;
    // Setup the workspace for the method
    int limit = 10000;
    int lenw = 100000;
    int last = 0;
    int * iwork = new int[limit];
    double * work = new double[lenw];
    // Evaluate the integral
    dqags_(fF2Photon, params, &a, &b, &epsabs, &epsrel, &izn, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    izn = std::pow(Q2, J) * izn / alphaQED(Q2);
    // Free workspace from memory
    delete[] iwork;
    delete[] work;
    return izn;
}

double F2Photon::IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs)
{
    /*
        Computes the IzNBar integral that appears in the 
        holographic computation of F2Photon
        kin - std::vector<double> with just one element: W
        reg - Reggeon object from wich the reggeon index can be accessed.
        gs - std::vector<double> that contains the constant quantities in our problem.
        returns x^(1-j_n) g_n / (4 pi alphaE(0)) associated with reggeon n
    */
   const double alpha0 = 0.0072973525693;
   const double x = kin[0];
   const double J = reg.getJ();
   const int reg_index = reg.getIndex();
   double iznbar = std::pow(x, 1-J) * gs[reg_index-1] / (4 * M_PI * alpha0);
   return iznbar;
}

std::vector<kinStruct>  F2Photon::getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec)
{
    // Create unique list of Q2s
    std::vector<double> Q2s = points[0];
    sort(Q2s.begin(), Q2s.end());
    Q2s.erase(unique(Q2s.begin(), Q2s.end()), Q2s.end());
    const int n_Q2s = Q2s.size();
    // Get the Reggeons in spec
    const std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    // Compute the modes first if not computed previously
    if(modes.size() == 0) this->computeU1NNModes();
    // Create vector of kinStructs
    std::vector< kinStruct > ans(n_Q2s);
    std::vector<double> kinematics, izns(n_reggeons,0);
    for(int i = 0; i < n_Q2s; i++)
    {
        kinematics = {Q2s[i]};
        for(int reg_idx = 0; reg_idx < n_reggeons; reg_idx++) izns[reg_idx] = IzN(kinematics, reggeons[reg_idx]);
        kinStruct izs(Q2s[i], izns);
        ans[i] = izs;
    }
    return ans ;
}

std::vector<kinStruct>  F2Photon::getIzsBar(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec,
                                          const std::vector<double>  &gs)
{
    /*
        Computes all the necessary IzNBar.
        points: std::vector of vectors. Each std::vector contains a pair (Q2, x).
        spec: object which contains a std::vector with different kernels.
        gs: std::vector with the gn corresponding to each reggeon present in spec.
        To get the correct gs for each point we initiate the variable gs_index.
        Then it iterates over the kernels and then over the reggeons of each kernel.
        In each iteration gs_index is increased by one.
        This guarantees that the correct gn is obtained.        
    */
    // Create list of unique xs
    std::vector<double> xs = points[1];
    std::sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());
    const int n_xs = xs.size();
    std::vector< kinStruct > ans(n_xs) ;
    // Get kernels in spec. For F2 there is only one value of t, i.e. 0
    std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    std::vector<double> kinematics, iznbars(n_reggeons,0);
    for(int i = 0; i < n_xs; i++)
    {
        kinematics = {xs[i]};
        for(int reg_index = 0; reg_index < n_reggeons; reg_index++) iznbars[reg_index] = IzNBar(kinematics, reggeons[reg_index], gs);
        kinStruct izbars(xs[i], iznbars);
        ans[i] = izbars;
    }
    return ans ;
}

std::vector<double>  F2Photon::predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar, const std::vector< std::vector<double> >  &points, const bool savePredictions)
{
    // Check that Izs, IzsBar and points have the same size
    if (points.size() == 0) throw std::runtime_error("points has length 0. Aborting F2Photon::predict.");
    std::vector<double> ans(points[0].size()) ;
    kinStruct iznStruct, iznbarStruct;
    std::vector<double> izn, iznbar;
    for(int i = 0; i < points[0].size(); i++)
    {
        iznStruct = binary_search<kinStruct>(Izs, kinStruct(points[0][i],{}));
        iznbarStruct = binary_search<kinStruct>(IzsBar, kinStruct(points[1][i],{}));
        izn = iznStruct.izns; iznbar = iznbarStruct.izns;
        ans[i] = sum(izn * iznbar);
    }
    if(savePredictions)
    {
        std::ofstream myfile;
        std::string file_path;
        std::cout << "Please introduce the path to save the predictions of F2Photon" << std::endl;
        std::cin >> file_path;
        myfile.open(file_path);
        myfile << "Q2\tx\tPred" << std::endl;
        for(int i = 0; i < points[0].size(); i++)
        {
            myfile << points[0][i] << '\t' << points[1][i] << '\t' << ans[i] << std::endl;
        }
    }
    return ans ;
}

F2Photon& F2Photon::operator=(const F2Photon &rhs)
{
    if (this == &rhs) return *this;
    copy(rhs);
    return *this;
}

F2Photon::~F2Photon()
{
    // Class destructor
}