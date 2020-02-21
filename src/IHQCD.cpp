#include "IHQCD.h"
#include <boost/numeric/odeint.hpp>
#include "methods/vectorOperators.hpp"
#include "methods/interpolation/Spline_Interp.hpp"

const double IHQCD::b0 = 4.2;
const double IHQCD::b1 = 51 * std::pow(IHQCD::b0, 2) / 121.0 ;
const double IHQCD::alpha = 2;
const double IHQCD::aa = (3.0/8) * (IHQCD::alpha - 1) / IHQCD::alpha ;
const double IHQCD::LAdS = 1;

// Constructor with arguments
IHQCD::IHQCD(const double a0, const double l0, const double z, const double delta):
    A0(a0), lambda0(l0), zmax(z), h(delta),
    zs({}), As({}),
    dAs({}), d2A({}), d3A({}),
    lambdas({}), dlambdas({}), d2lambda({}), d3lambda({}),
    Astrings({}), dAstrings({}), d2Astring({}), d3Astring({}),
    Phis({}), dPhis({}), d2Phi({}), d3Phi({}),
    u0({}), u2({}),
    aF({}), bF({}), cF({}), l1_2({}),
    e2As({}), e2A({})
    {
        solve();
    }

IHQCD::IHQCD(const IHQCD &ihqcd):
    A0(ihqcd.A0), lambda0(ihqcd.lambda0),
    zmax(ihqcd.zmax), h(ihqcd.h),
    zs(ihqcd.zs), As(ihqcd.As),
    dAs(ihqcd.dAs), d2A(ihqcd.d2A), d3A(ihqcd.d3A),
    lambdas(ihqcd.lambdas), dlambdas(ihqcd.dlambdas), d2lambda(ihqcd.d2lambda), d3lambda(ihqcd.d3lambda),
    Astrings(ihqcd.Astrings), dAstrings(ihqcd.dAstrings), d2Astring(ihqcd.d2Astring), d3Astring(ihqcd.d3Astring),
    Phis(ihqcd.Phis), dPhis(ihqcd.dPhis), d2Phi(ihqcd.d2Phi), d3Phi(ihqcd.d3Phi),
    u0(ihqcd.u0), u2(ihqcd.u2),
    aF(ihqcd.aF), bF(ihqcd.bF), cF(ihqcd.cF), l1_2(ihqcd.l1_2),
    e2As(ihqcd.e2As), e2A(ihqcd.e2A)
    {}

void IHQCD::copy(const IHQCD &ihqcd)
{
    A0 = ihqcd.A0; lambda0 = ihqcd.lambda0;
    zmax = ihqcd.zmax; h = ihqcd.h;
    zs = ihqcd.zs;
    As = ihqcd.As; dAs = ihqcd.dAs; d2A = ihqcd.d2A; d3A = ihqcd.d3A;
    lambdas = ihqcd.lambdas; dlambdas = ihqcd.dlambdas; d2lambda = ihqcd.d2lambda; d3lambda = ihqcd.d3lambda;
    Astrings = ihqcd.Astrings; dAstrings = ihqcd.dAstrings; d2Astring = ihqcd.d2Astring; d3Astring = ihqcd.d3Astring;
    Phis = ihqcd.Phis; dPhis = ihqcd.dPhis; d2Phi = ihqcd.d2Phi; d3Phi = ihqcd.d3Phi;
    u0 = ihqcd.u0; u2 = ihqcd.u2;
    aF = ihqcd.aF; bF = ihqcd.bF; cF = ihqcd.cF; l1_2 = ihqcd.l1_2;
    e2As = ihqcd.e2As; e2A = ihqcd.e2A;
}

void IHQCD::finalizeBackground()
{
    zs = std::vector<double>(zs.begin() + 5, zs.end());
    As = std::vector<double>(As.begin() + 5, As.end());
    dAs = std::vector<double>(dAs.begin() + 5, dAs.end());
    d2A = std::vector<double>(d2A.begin() + 5, d2A.end());
    d3A = std::vector<double>(d3A.begin() + 5, d3A.end());
    lambdas = std::vector<double>(lambdas.begin() + 5, lambdas.end());
    dlambdas = std::vector<double>(dlambdas.begin() + 5, dlambdas.end());
    d2lambda = std::vector<double>(d2lambda.begin() + 5, d2lambda.end());
    d3lambda = std::vector<double>(d3lambda.begin() + 5, d3lambda.end());
    Astrings = std::vector<double>(Astrings.begin() + 5, Astrings.end());
    dAstrings = std::vector<double>(dAstrings.begin() + 5, dAstrings.end());
    d2Astring = std::vector<double>(d2Astring.begin() + 5, d2Astring.end());
    d3Astring = std::vector<double>(d3Astring.begin() + 5, d3Astring.end());
    Phis = std::vector<double>(Phis.begin() + 5, Phis.end());
    dPhis = std::vector<double>(dPhis.begin() + 5, dPhis.end());
    d2Phi = std::vector<double>(d2Phi.begin() + 5, d2Phi.end());
    d3Phi = std::vector<double>(d3Phi.begin() + 5, d3Phi.end());  
    u0 = std::vector<double>(u0.begin() + 5, u0.end());
    u2 = std::vector<double>(u2.begin() + 5, u2.end());
    aF = std::vector<double>(aF.begin() + 5, aF.end());
    bF = std::vector<double>(bF.begin() + 5, bF.end());
    cF = std::vector<double>(cF.begin() + 5, cF.end());
    l1_2 = std::vector<double>(l1_2.begin() + 5, l1_2.end());
    e2As = std::vector<double>(e2As.begin() + 5, e2As.end());
    e2A = std::vector<double>(e2A.begin() + 5, e2A.end());
}

double IHQCD::X(const double l)
{
    // Returns the superpotential X(\lambda)
    double ans = - b0 * l / (3 + 2 * b0 * l) - (2 * b0 * b0 + 3 * b1) * l * l / (9 * (1 + l * l ) * (1 + (2 * b0 * b0 + 3 * b1) * std::log(1 + l * l) / (18 * aa))) ;
    return ans;
}

double IHQCD::W(const double l)
{
    // Returns the superpotential W(\l)
    double ans = (9./4) * std::pow(1 + (2./3) * b0 * l, 2.0/3) * std::pow(1 + (2 * b0 * b0 + 3 * b1) * std::log(1 + l * l) / (18. * aa), 4 * aa / 3.0);
    return ans;
}

double IHQCD::dAEOM(const double A, const double l)
{
    // Computes dAdz given A and lambda using the EOMs
    const double WW = W(l);
    return - (4.0/9) * WW * std::exp(A) / LAdS ;
}

double IHQCD::dlambdaEOM(double A, double l)
{
    // Computes dAdz given A and lambda using the EOMs
    const double XX = X(l);
    const double WW = W(l);
    return - (4.0/3) * XX * WW * l * std::exp(A) / LAdS ;
}

void IHQCD::eom(const state &Y , state &dYdz , const double z)
{
    // Y = (A, lambda, dA, dlambda)
    const double A = Y[0] , lambda = Y[1] ;
    const double XX = X(lambda);
    const double WW = W(lambda);
    dYdz[0] = - (4.0/9) * WW * std::exp(A) / LAdS ;
    dYdz[1] = - (4.0/3) * XX * WW * lambda * std::exp(A) / LAdS ;
}

void IHQCD::observer(const state &Y, const double z)
{
    // Y = [A, lambda]
    const double l = Y[1];
    const double phi = std::log(l);
    const double da =  dAEOM(Y[0], l);
    const double dl = dlambdaEOM(Y[0], l);
    const double dphi = dl / l;
    const double as = Y[0] + (2.0 / 3) * phi ;
    const double das = da + (2.0/3) * dphi;
    zs.push_back(z);
    As.push_back(Y[0]);
    lambdas.push_back(l);
    dAs.push_back(da);
    dlambdas.push_back(dl);
    Phis.push_back(phi);
    dPhis.push_back(dphi);
    Astrings.push_back(as);
    dAstrings.push_back(das);
    
    // cF, l1_2, e2As, e2A
    cF.push_back(dphi * dphi);
    l1_2.push_back(std::sqrt(l)); e2As.push_back(std::exp(2 * as)); e2A.push_back(std::exp(2 * Y[0]));
}

void IHQCD::solve()
{
    /* 
        Solves the Gursoy, Kiritsis and Nitti model in the Einstein frame
        Finds A and lambda as a function of z.
    */
   std::cout << "Solving IHQCD for A0 = " << A0 << ", lambda0 = " << lambda0 << ", zmax = " << zmax << ", h = " << h << std::endl;
    // zmin = exp(-A0)
    const double z = std::exp(-A0);
    // Boundary conditions
    state XYM = {A0, lambda0};
    auto eomfun = [this] (const state &Y, state &dYdA, double A) {this->eom(Y, dYdA, A);};
    auto obsfun = [this] (const state &Y, double A) {this->observer(Y, A);};
    // Integration of the EOM
    boost::numeric::odeint::dense_output_runge_kutta< boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<state > > > stepper;
    integrate_const(stepper, eomfun, XYM, z, zmax, h, obsfun);
    // Now we compute d2A and d3A
    Spline_Interp<double>  dAProfile(zs, dAs);
    for( int i = 0; i < zs.size(); i++)
    {
        d2A.push_back(dAProfile.der1(zs[i]));
        d3A.push_back(dAProfile.der2(zs[i]));
    }

    // Now we compute d2lamba and d3lambda
    Spline_Interp<double>  dlambdaProfile(zs, dlambdas);
    for( int i = 0; i < zs.size(); i++)
    {
        d2lambda.push_back(dlambdaProfile.der1(zs[i]));
        d3lambda.push_back(dlambdaProfile.der2(zs[i]));
    }

    // Now we compute d2phi and d3phi
    Spline_Interp<double>  dphiProfile(zs, dPhis);
    for( int i = 0; i < zs.size(); i++)
    {
        d2Phi.push_back(dphiProfile.der1(zs[i]));
        d3Phi.push_back(dphiProfile.der2(zs[i]));
    }

    // Computing d2Astring and d3Astring
    Spline_Interp<double>  dAstringProfile(zs, dAstrings);
    for( int i = 0; i < zs.size(); i++)
    {
        d2Astring.push_back(dAstringProfile.der1(zs[i]));
        d3Astring.push_back(dAstringProfile.der2(zs[i]));
    }
    // Scalar potential
    u0 = (9.0/4) * dAs * dAs - 1.5 * d2A + 2.0 * d2A * d2A / (dAs * dAs) + 3.0 * dAs * d2Phi / dPhis - 2.0 * d2A * d2Phi / (dAs * dPhis) - d3A / dAs + d3Phi / dPhis ;
    // Tensor potential
    u2 = (15.0/4) * dAs * dAs - (2.0/3) * dPhis * dPhis;
    
    // aF, bF
    aF = d2Phi;
    bF = d2Astring - dAstrings * dAstrings;
    finalizeBackground();
}

std::vector<double> IHQCD::getZ()
{
    return zs;
}

std::vector<double> IHQCD::getA()
{
    return As;
}

std::vector<double> IHQCD::getdA()
{
    return dAs;
}

std::vector<double> IHQCD::getd2A()
{
    return d2A;
}

std::vector<double> IHQCD::getd3A()
{
    return d3A;
}

std::vector<double> IHQCD::getlambda()
{
    return lambdas;
}

std::vector<double> IHQCD::getdlambda()
{
    return dlambdas;
}

std::vector<double> IHQCD::getd2lambda()
{
    return d2lambda;
}

std::vector<double> IHQCD::getd3lambda()
{
    return d3lambda;
}

std::vector<double> IHQCD::getAstring()
{
    return Astrings;
}

std::vector<double> IHQCD::getdAstring()
{
    return dAstrings;
}

std::vector<double> IHQCD::getd2Astring()
{
    return d2Astring;
}

std::vector<double> IHQCD::getd3Astring()
{
    return d3Astring;
}

std::vector<double> IHQCD::getPhi()
{
    return Phis;
}
std::vector<double> IHQCD::getdPhi()
{
    return dPhis;
}

std::vector<double> IHQCD::getd2Phi()
{
    return d2Phi;
}

std::vector<double> IHQCD::getd3Phi()
{
    return d3Phi;
}

std::vector<double> IHQCD::getU0()
{
    return u0;
}

std::vector<double> IHQCD::getU2()
{
    return u2;
}

std::vector<double> IHQCD::getaF()
{
    return aF;
}

std::vector<double> IHQCD::getbF()
{
    return bF;
}

std::vector<double> IHQCD::getcF()
{
    return cF;
}

std::vector<double> IHQCD::getl1_2()
{
    return l1_2;
}

std::vector<double> IHQCD::getE2As()
{
    return e2As;
}

std::vector<double> IHQCD::getE2A()
{
    return e2A;
}

IHQCD& IHQCD::operator= (const IHQCD &rhs)
{
    if (this == &rhs) return *this;
    copy(rhs);
    return *this;
}

IHQCD::~IHQCD(){}

IHQCD& ihqcd()
{
    static IHQCD bck(5, 0.0337462, 10.0, 0.002);
    return bck;
}

IHQCD& ihqcdU1NNMode()
{
    static IHQCD bck(5, 0.0337462, 20.0, 0.001);
    return bck;
}