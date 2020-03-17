#include <iostream>
#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>

#include "IHQCD.h"
#include "SigmaGammaGamma.h"
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

class SigmaIzNIntegrand
{
    private:
        Poly_Interp<double> f1, f2;
    public:
        SigmaIzNIntegrand(const Poly_Interp<double> &func1, const Poly_Interp<double> &func2);
        double operator()(const double x);
};

SigmaIzNIntegrand::SigmaIzNIntegrand(const Poly_Interp<double> &func1, const Poly_Interp<double> &func2): f1(func1), f2(func2) {}

double SigmaIzNIntegrand::operator()(const double x) {return f1.interp(x) * f2.interp(x);}

double fSigma(double * x, void * params)
{
    return ((SigmaIzNIntegrand *) params)->operator()(*x);
}

void SigmaGammaGamma::loadData(std::string file_path)
{
    /*
        This function loads the data of sigma(gamma gamma -> hadrons) contained in the path file_path
        The data is assumed to be in 6 columns:
        W(GeV), WErrorp, WErrorm, sigma(np), Errorp and Errorm
        WErrorp is the positive error uncertainty in W
        WErrorm is the negative error uncertainty in W
        sigma(nb) is the cross section in nb
        Errorp is the plus error in sigma
        Errorm is the minus error in sigma 
    */
   // Check that the file path is not empty
   if (file_path.size() == 0)
   {
       std::cout << "No file path given. Exiting SigmaGammaGamma::loadData" << std::endl;
       return;
   }
   std::ifstream file;
   file.open(file_path);
   if(file.fail()) std::runtime_error("Error: File with sigma(gamma gamma -> hadrons) data not opened."); 
   std::string line ;
   getline(file, line) ;
   std::cout << "Loading sigma(gamma gamma -> hadrons) data" << std::endl;
   std::vector<double> Ws, sigmas, sigmaErrs;
   std::vector<std::string> result;
   while(getline(file, line))
    {
        boost::split(result, line, boost::is_any_of("\t") ) ;
        Ws.push_back(stod(result[0])) ;
        sigmas.push_back(stod(result[3]));
        sigmaErrs.push_back(std::max(stod(result[4]), stod(result[5])));
    }
    // Check that all std::vector containers have the same side
    if( Ws.size() != sigmas.size() || Ws.size() != sigmaErrs.size() || sigmas.size() != sigmaErrs.size())
    {
        throw std::runtime_error("Ws, sigma and sigmaErrs vector containers in SigmaGammaGamma don't have the same size"); 
    }
    this->setDataPts({Ws, sigmas, sigmaErrs});
}

SigmaGammaGamma::SigmaGammaGamma(std::string file_path): PhotonScattering()
{
    loadData(file_path);
}

std::vector<double> SigmaGammaGamma::expVal()
{
    // Returns a std::vector container with the values of sigma
    return this->getDataPts()[1];
}

std::vector<double>  SigmaGammaGamma::expErr()
{
    // Returns a std::vector container with the values of sigmaErrs
    return this->getDataPts()[2];
}

std::vector<std::vector<double> >  SigmaGammaGamma::expKinematics()
{
    /*
        Returns a std::vector<std::vector<double> > with the values of Ws.
        This format is just to make it uniform with processed that depend
        in more than one kinematical variable
    */
   std::vector<double> Ws = this->getDataPts()[0];
   return {Ws};
}

double SigmaGammaGamma::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    /*
        Computes the IzN integral that appears in the 
        holographic computation of sigma(gamma gamma -> hadrons)
        kin - std::vector<double> with just one element: W
        reg - Reggeon object from wich j_n and \psi_n can be accessed.
        returns izn = \int dz e^{-(j_n - 1.5) A_s} \psi_n
    */
    // Get J and the wavefunction wf from Reggeon object reg
    const double J = reg.getJ();
    std::vector<std::vector<double> > wf = reg.getWf();

    // Compute e^{-(j_n - 1.5) A_s}
    std::vector<double> z = ihqcd().getZ();
    std::vector<double> fact1 = exp((1.5-J) * ihqcd().getAstring());
    Poly_Interp<double> f1(z, fact1, 4);

    // Compute \psi_n
    Poly_Interp<double> f2(wf[0], wf[1], 4);

    // Define the integrand object
    SigmaIzNIntegrand integrand(f1, f2);
    void * params = &integrand;
    // Compute the integral
    double izn = 0.0, abserr = 0.0;
    double a = z[0], b = z.back();
    // Absolute and relative tolerances desidered
    double epsabs = 1e-9, epsrel = 1e-9;
    // Output variables 
    int neval = 0, ier;
    // Setup the workspace for the method
    int limit = 10000;
    int lenw = 100000;
    int last = 0;
    int * iwork = new int[limit];
    double * work = new double[lenw];
    // Evaluate the integral
    dqags_(fSigma, params, &a, &b, &epsabs, &epsrel, &izn, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    // Free workspace from memory
    delete[] iwork;
    delete[] work;
    return izn;
}

double SigmaGammaGamma::IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs)
{
    /*
        Computes the IzNBar integral that appears in the 
        holographic computation of sigma(gamma gamma -> hadrons)
        kin - std::vector<double> with just one element: W
        reg - Reggeon object from wich the reggeon index can be accessed.
        gs - std::vector<double> that contains the constant quantities in our problem.
        returns s^(j_n -1 ) g_n associated with reggeon n
    */
   // Compute s
   double s = std::pow(kin[0],2), J = reg.getJ();
   const int reg_index = reg.getIndex();
   double iznbar = std::pow(s, J - 1) * gs[reg_index-1];
   return iznbar;
}

std::vector<kinStruct>  SigmaGammaGamma::getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec)
{
    /*
        Computes all the Izs relevant to sigma(gamma gamma -> hadrons)
    */
    // Get the Reggeons in spec
    const std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    // Create vector of kinStructs
    std::vector< kinStruct > ans(1);
    std::vector<double> izns(n_reggeons,0);
    for(int reg_idx = 0; reg_idx < n_reggeons; reg_idx++) izns[reg_idx] = IzN({}, reggeons[reg_idx]);
    kinStruct izs(0.0, izns);
    ans[0] = izs;
    return ans ;
}

std::vector<kinStruct> SigmaGammaGamma::getIzsBar(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec,
                                          const std::vector<double>  &gs)
{
    /*
        Computes all the IzNbars relevant to sigma(gamma gamma -> hadrons)
    */
    // Create list of unique Ws
    std::vector<double> Ws = points[0];
    std::sort(Ws.begin(), Ws.end());
    Ws.erase(unique(Ws.begin(), Ws.end()), Ws.end());
    const int n_Ws = Ws.size();
    std::vector<kinStruct> ans(n_Ws) ;
    // Get kernels in spec. For sigma(gamma gamma -> hadrons) there is only one value of t, i.e. 0
    std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    std::vector<double> kinematics, iznbars(n_reggeons,0);
    for(int i = 0; i < n_Ws; i++)
    {
        kinematics = {Ws[i]};
        for(int reg_index = 0; reg_index < n_reggeons; reg_index++) iznbars[reg_index] = IzNBar(kinematics, reggeons[reg_index], gs);
        kinStruct izbars(Ws[i], iznbars);
        ans[i] = izbars;
    }
    return ans ;
}

std::vector<double>  SigmaGammaGamma::predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                     const std::vector< std::vector<double> >  &points, const bool savePredictions)
{
    // Check that Izs, IzsBar and points have the same size
    if (points.size() == 0) throw std::runtime_error("points has length 0. Aborting SigmaGammaGamma::predict.");
    std::vector<double> ans(points[0].size()) ;
    kinStruct iznStruct, iznbarStruct;
    std::vector<double> izn, iznbar;
    for(int i = 0; i < points[0].size(); i++)
    {
        iznStruct = binary_search<kinStruct>(Izs, kinStruct({0.0},{}));
        iznbarStruct = binary_search<kinStruct>(IzsBar, kinStruct(points[0][i],{}));
        izn = iznStruct.izns; iznbar = iznbarStruct.izns;
        ans[i] = sum(izn * iznbar);
    }
    if(savePredictions)
    {
        std::ofstream myfile;
        std::string file_path;
        std::cout << "Please introduce the path to save the predictions of sigma(gamma gamma -> hadrons)" << std::endl;
        std::cin >> file_path;
        myfile.open(file_path);
        myfile << "W\tPred" << std::endl;
        for(int i = 0; i < points[0].size(); i++) myfile << points[0][i] << '\t' << ans[i] << std::endl;
    }
    return ans ;
}

SigmaGammaGamma::~SigmaGammaGamma()
{
    // Class destructor
}