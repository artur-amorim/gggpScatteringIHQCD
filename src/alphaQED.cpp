#include "alphaQED.h"

#include <cmath>

// Definition of the fine structure constant, i.e alphaQED(0)
const double alpha0 = 0.0072973525693;
// Definition of the electron mass in GeV
const double me = 0.0005109989500;
// Definition of the muon mass in GeV
const double mmu = 0.1056583755;
// Definition of the tau mass in GeV
const double mtau = 1.77686;
// Definition of the Z boson mass in GeV
const double mZ = 91.1876;
// Definition of Riemann zeta(2)
const double ZETA_2 = 1.64493406684823;
// Definition of Riemann zeta(3)
const double ZETA_3 = 1.20205690315959;
// Definition of Riemann zeta(5)
const double ZETA_5 = 1.03692775514337;

double Pi0(const double s, const double m1)
{
    /* 
        One loop contribution of a lepton of mass m1 to alphaQED.
        Formula taken from
        "Leptonic contribution to the effective electromagnetic coupling constant up to three loops" by M. Steinhauser, arXiv hep-ph 9803313
    */
   double ans = 20./9 - (4.0/3) * std::log(s/(m1*m1)) + 8 * m1 * m1 / s ;
   return ans;
}

double Pi1(const double s, const double m1)
{
    /* 
        Two loop contribution of a lepton of mass m1 to alphaQED.
        Formula taken from
        "Leptonic contribution to the effective electromagnetic coupling constant up to three loops" by M. Steinhauser, arXiv hep-ph 9803313
    */
   double Lqm1 = std::log(s/(m1*m1));
   double ans = 5./6 - 4 * ZETA_3 - Lqm1 - 12 * Lqm1 * m1 * m1 / s ;
   return ans;
}

double Pi2A(const double s, const double m1)
{
    /* 
        Three loop quenched contribution of a lepton of mass m1 to alphaQED.
        Formula taken from
        "Leptonic contribution to the effective electromagnetic coupling constant up to three loops" by M. Steinhauser, arXiv hep-ph 9803313
    */
   double ans = -(121./48) + (-5 + 8 * std::log(2)) * ZETA_2 - (99./16) * ZETA_3 + 10 * ZETA_5 + std::log(s/(m1*m1)) / 8.0 ;
   return ans;
}

double Pi2l(const double s, const double m1, const double m2)
{
    /* 
        Three loop Pi2l contribution of a bubble diagram.
        For more details see
        "Leptonic contribution to the effective electromagnetic coupling constant up to three loops" by M. Steinhauser, arXiv hep-ph 9803313
    */
   const double Lqm1 = std::log(s / (m1 * m1)), Lqm2 = std::log(s / (m2 * m2));
   double ans = - 116./27 + (4.0/3) * ZETA_2 + (38./9) * ZETA_3 + (14./9) * Lqm1 + ((5./18) - (4./3) * ZETA_3) * Lqm2 + (1./6) * Lqm1 * Lqm1 - Lqm1 * Lqm2 / 3.0;
   return ans;
}

double Pi2F(const double s, const double m1)
{
    /* 
        Three loop Pi2f contribution of a bubble diagram.
        For more details see
        "Leptonic contribution to the effective electromagnetic coupling constant up to three loops" by M. Steinhauser, arXiv hep-ph 9803313
    */
   const double Lqm1 = std::log(s / (m1 * m1));
   double ans = - 307./216 - (8.0/3) * ZETA_2 + (545./144) * ZETA_3 + (11./6 - (4./3) * ZETA_3) * Lqm1 - Lqm1 * Lqm1 / 6.0;
   return ans;
}

double Pi2h(const double s, const double m2)
{
    /* 
        Three loop Pi2h contribution of a bubble diagram.
        For more details see
        "Leptonic contribution to the effective electromagnetic coupling constant up to three loops" by M. Steinhauser, arXiv hep-ph 9803313
    */
   const double Lqm2 = std::log(s / (m2 * m2));
   double ans = - 37./6 + 38 * ZETA_3 / 9 + (11./6 - 4 * ZETA_3 / 3) * Lqm2 - Lqm2 * Lqm2 / 6 ;
   return ans;
}

double DeltaAlphaLep(const double s)
{
    /* Evaluates the contribution of the leptons to the photon propagator
        After renormalization this contribution induces a running of alphaQED
        This function takes s or t in GeV^2.
        We follow the computation outlined in
        "Leptonic contribution to the effective electromagnetic coupling constant up to three loops" by M. Steinhauser, arXiv hep-ph 9803313
    */
   // Compute one loop corrections
   double Pi0_electron = Pi0(s, me);
   double Pi0_muon = Pi0(s, mmu);
   double Pi0_tau = Pi0(s, mtau);
   double Delta_Alpha_0 = - alpha0 * (Pi0_electron + Pi0_muon + Pi0_tau) / (4 * M_PI);

   // Compute two loop corrections
   double Pi1_electron = Pi1(s, me);
   double Pi1_muon = Pi1(s, mmu);
   double Pi1_tau = Pi1(s, mtau);
   double Delta_Alpha_1 = - (alpha0 * alpha0 / (4 * M_PI * M_PI)) * (Pi1_electron + Pi1_muon + Pi1_tau) ;

   // Compute three loop corrections
   // Quenched contribution
   double Pi2A_electron = Pi2A(s, me);
   double Pi2A_muon = Pi2A(s, mmu);
   double Pi2A_tau = Pi2A(s, mtau);
   double Delta_Alpha_2 = Pi2A_electron + Pi2A_muon + Pi2A_tau;
   // Add Pi2l contributions
   Delta_Alpha_2 += Pi2l(s, mmu, me) + Pi2l(s, mtau, me) + Pi2l(s,mtau, mmu);
   // Add Pi2F constributions
   Delta_Alpha_2 += Pi2F(s, me) + Pi2F(s, mmu) + Pi2F(s, mtau);
   // Add Pi2h contributions
   Delta_Alpha_2 += Pi2h(s, mmu) + 2 * Pi2h(s, mtau);
   Delta_Alpha_2 = - 0.25 * std::pow(alpha0/M_PI,3) * Delta_Alpha_2;
   // Return total DeltaAlphaLeptons
   double ans =  Delta_Alpha_0 + Delta_Alpha_1 + Delta_Alpha_2;
   return ans;
}

double DeltaAlphaHad(const double s)
{
    /* Evaluates the hadronic contribution to the photon propagator
        After renormalization this contribution induces a running of alphaQED
        This function takes s or t in GeV^2.
        This is a parametrization proposed in the paper
        "Update of the hadronic contribution to the QED vacuum polarization"
    */
   double A, B, C, sqrtS = std::sqrt(s);

   // Values of A, B and C are present in the table 3 of the paper mentioned above
   if (sqrtS <= 0.7)
   {
       A = 0.0;
       B = 0.0023092;
       C = 3.9925370;
   }
   else if ( sqrtS > 0.7 && sqrtS <= 2.0)
   {
       A = 0.0;
       B = 0.0022333;
       C = 4.2191779;
   }
   else if (sqrtS > 2.0 && sqrtS <= 4.0)
   {
       A = 0.0;
       B = 0.0024402;
       C = 3.2496684;
   }
   else if (sqrtS > 4.0 && sqrtS <= 10)
   {
       A = 0.0;
       B = 0.0027340;
       C = 2.0995092;
   }
   else if (sqrtS > 10.0 && sqrtS <= mZ)
   {
       A = 0.0010485;
       B = 0.0029431;
       C = 1.0;
   }
   else if (sqrtS > mZ && sqrtS <= 10000.0)
   {
       A = 0.0012234;
       B = 0.0029237;
       C = 1.0;
   }
   else
   {
       A = 0.0016894;
       B = 0.0028984;
       C = 1.0;
   }
   
    double DeltaAlpha = A + B * std::log(1 + C * std::fabs(s));
    return DeltaAlpha;
}

double alphaQED(const double s)
{
    /*  Computes alphaQED = alpha / (1 - DeltaAlphaLep - DeltaAlphaHad)
        s is given in GeV^2
    */
   double delta_alpha_leptons = DeltaAlphaLep(s);
   double delta_alpha_quarks = DeltaAlphaHad(s);
   double alpha = alpha0 / (1 - delta_alpha_leptons - delta_alpha_quarks);
   return alpha;
}

double InvalphaQED(const double s)
{
    /*  Computes alphaQED = alpha / (1 - DeltaAlphaLep - DeltaAlphaHad)
        s is given in GeV^2
    */
   double delta_alpha_leptons = DeltaAlphaLep(s);
   double delta_alpha_quarks = DeltaAlphaHad(s);
   double inv_alpha = (1 - delta_alpha_leptons - delta_alpha_quarks) / alpha0;
   return inv_alpha;
}