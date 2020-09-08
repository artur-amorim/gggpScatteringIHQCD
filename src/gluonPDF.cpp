#include <fstream>
#include <cmath>
#include "gluonPDF.h"
#include "IHQCD.h"
#include "HardPomeron.h"
#include "methods/search.hpp"
#include "LHAPDF/LHAPDF.h"
#include "schrodinger/schrodinger.h"
#include "methods/vectorOperators.hpp"

GluonPDF::GluonPDF(): 
    f2(F2("")), fl(FL("")), pomeron(HQCDP())
    {
        HardPomeron * gluon = new HardPomeron(4, {});
        pomeron.addKernel(*gluon);
        // Add f2 and fl to HQCDP object
        pomeron.addProcessObservable(f2);
        pomeron.addProcessObservable(fl);
    }

GluonPDF::GluonPDF(const GluonPDF &rhs): f2(rhs.f2), fl(rhs.fl), pomeron(rhs.pomeron) {}

double alphaSQCD(const double Q2)
{
    /*
        Returns the alphaS of QCD using AlphaS_Ipol object from LHAPDF
        and the given values of alpha_S and corresponding Q in the PDF set
        CT18ZNNLO
    */
    LHAPDF::AlphaS_Ipol as_ipol_CT18;

    std::vector<double> AlphaS_Qs = {1.39500E+00, 1.39875E+00, 1.59202E+00, 1.82174E+00, 2.09834E+00, 2.43364E+00, 2.84296E+00, 3.34637E+00, 3.97032E+00,
                 4.75000E+00, 5.76715E+00, 7.07072E+00, 8.75819E+00, 1.09657E+01, 1.38856E+01, 1.77929E+01, 2.30855E+01, 3.03471E+01,
                 4.04448E+01, 5.46864E+01, 7.50724E+01, 1.04712E+02, 1.48517E+02, 2.14380E+02, 3.15212E+02, 4.72537E+02, 7.22946E+02,
                 1.12995E+03, 1.80616E+03, 2.95593E+03, 4.95886E+03, 8.53814E+03, 1.51079E+04, 2.75107E+04, 5.16275E+04, 1.00000E+05};
    std::vector<double> AlphaS_Vals = {3.62194E-01, 3.61558E-01, 3.37518E-01, 3.14861E-01, 2.94393E-01, 2.75800E-01, 2.58830E-01, 2.43276E-01, 2.28967E-01,
                   2.15760E-01, 2.04638E-01, 1.93923E-01, 1.83856E-01, 1.74388E-01, 1.65473E-01, 1.57071E-01, 1.49144E-01, 1.41660E-01,
                   1.34588E-01, 1.27902E-01, 1.21575E-01, 1.15586E-01, 1.09913E-01, 1.04537E-01, 9.94398E-02, 9.46049E-02, 9.00171E-02,
                   8.56622E-02, 8.15271E-02, 7.75995E-02, 7.38679E-02, 7.03216E-02, 6.69506E-02, 6.37456E-02, 6.06978E-02, 5.77990E-02};

    as_ipol_CT18.setQValues(AlphaS_Qs);
    as_ipol_CT18.setAlphaSValues(AlphaS_Vals);

    double alphaS = as_ipol_CT18.alphasQ2(Q2);
    return alphaS;
}

double alphaSIHQCD(const double E)
{
    /*
        Returns the alphaS prediction of IHQCD. In IHQCD is given by alphaS = c0 * lambda / 12 PI. 
    */
    const double b0 = 4.2;
    const double b0Bar = (2./3) * 11 / pow(4 * M_PI,2);
    const double c0 = b0 / b0Bar;

    const std::vector<double> z = ihqcd().getZ();
    const std::vector<double> As = ihqcd().getAstring();
    const std::vector<double> Phi = ihqcd().getPhi();
    const std::vector<double> lambda = ihqcd().getlambda();
    Poly_Interp<double> func(z, As - (2./3) * Phi, 4);
    Poly_Interp<double> lambFunc(z, lambda, 4);
    std::function<double(double)> f = [&func,E] (const double z) { return func.interp(z) - log(E);};
    double zRoot = zbrent(f, z[0], z.back(), 1e-9, false);
    double alphaS = c0 * lambFunc.interp(zRoot) / (12 * M_PI) ;
    return alphaS;
}

std::vector<double> GluonPDF::computePDF(const std::vector<std::vector<double> > & kinPts, const std::vector<double> pars)
{
    /*
     * Computes the PDF of the gluon for the kinematical points kinPts using kernel parameters and gn values in pars
     * Assumes kinPts = {Q2s, xs} and pars = {invls, a, b, c, d, g1, g2, g3, g4}
     * First it computes the spectrum after reading the values of invls, a, b, c and from pars
     * Next it computes the necessary points for F2 and FL as well their U1 NN modes
     * Then it computes the values of F2 and FL using the gns also present in pars
    */
    // Setup F2 and FL points
    f2.setDataPts(kinPts); fl.setDataPts(kinPts);
    // Compute the spectrum with invls, a, b, c, d give in pars
    pomeron.computeSpectrum({{pars[0], pars[1], pars[2], pars[3], pars[4]}});
    std::vector<Spectra> spec = pomeron.getSpectrum();
    // We are dealing with DIS observables so spec only has reggeons for t = 0
    std::vector<Reggeon> reggeons = spec[0].getReggeons();
    // gns needed to compute F2 and FL
    std::vector<double> gs = {pars[5], pars[6], pars[7], pars[8]};
    
    // Compute the modes of F2 and FL
    f2.computeU1NNModes();
    fl.computeU1NNModes();
    
    // Compute IzNBars. Can use f2 or fl object for that
    std::vector<kinStruct> xg_IzBars = f2.getIzsBar(kinPts, spec, gs);
    // Compute F2 IzN
    std::vector<kinStruct> f2Izs = f2.getIzs(kinPts, spec);
    // Compute FL IzN
    std::vector<kinStruct> flIzs = fl.getIzs(kinPts, spec);
    // Compute xg Izs
    std::vector<kinStruct> xg_Izs(f2Izs.size());
    const double Mc = 1.4;
    const double Mb = 4.75;
    const double muThc2 = Mc*Mc;
    const double muThb2 = Mb*Mb;

    // Now let's rescale the f2Izs and flIzs with the e2, alphaS and js factors
    double sum_ei2 = 0;
    double Q2 = 0;
    double alphaS = 0;
    double jn = 0;
    std::vector<double> izs(reggeons.size());
    for(int i = 0; i < f2Izs.size(); i++)
    {
        Q2 = f2Izs[i].kinVar;
        // u, d, s, c and b quarks are active
        if (Q2 >= muThb2) sum_ei2 = 11./9;
        // u, d, s, and c quarks are active
        else if ( Q2 >= muThc2 && Q2 < muThb2) sum_ei2 = 10./9;
        // only u, d, s are active
        else sum_ei2 = 2.0/3;
        alphaS = alphaSIHQCD(std::sqrt(Q2));
        for(int j = 0; j < reggeons.size(); j++)
        {
            jn = reggeons[j].getJ();
            izs[j] = f2Izs[i].izns[j] * (-4./3 - 2 * jn / 3) / sum_ei2;
            izs[j] += M_PI * flIzs[i].izns[j] * (1 + 1.5 * jn + 0.5 * jn * jn) / ( sum_ei2 * alphaS );
        }
        xg_Izs[i].kinVar = Q2;
        xg_Izs[i].izns = izs;
    }
    // Ok now it's time to compute xg(x, Q2)
    std::vector<double> xgs(kinPts[0].size(), 0.0);
    kinStruct iznStruct, iznbarStruct;
    std::vector<double> izn, iznbar;
    for(int i = 0; i < kinPts[0].size(); i++)
    {
        iznStruct = binary_search<kinStruct>(xg_Izs, kinStruct(kinPts[0][i],{}));
        iznbarStruct = binary_search<kinStruct>(xg_IzBars, kinStruct(kinPts[1][i],{}));
        izn = iznStruct.izns; iznbar = iznbarStruct.izns;
        xgs[i] = sum(izn * iznbar);
    }
    return xgs;
}

GluonPDF& GluonPDF::operator=(const GluonPDF &rhs)
{
    if (this == &rhs) return *this;
    f2 = rhs.f2;
    fl = rhs.fl;
    //pomeron = rhs.pomeron;
    return *this;
}

GluonPDF::~GluonPDF()
{
    // We need to delete the GluonKernel created in the constructor
    delete pomeron.getKernels()[0];
}