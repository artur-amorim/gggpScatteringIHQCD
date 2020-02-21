#ifndef IHQCD_H
#define IHQCD_H

#include <vector>

class IHQCD
{ 
    private:
        double A0, lambda0;             // Initial conditions used in solving the EOMs
        double zmax, h;                 // Maximum value of z coordinate and step used
        static const double b0, b1;
        static const double alpha;
        static const double aa;
        static const double LAdS;
        // Function X(\lambda)
        double X(const double l);
        // Function W(\lambda)
        double W(const double l);
        // dA / dz
        double dAEOM(const double A, const double l);
        // dlambda/dz
        double dlambdaEOM(double A, double l);
        // The operator() is used to define the EOMs
        typedef std::vector<double> state;
        void eom(const state &Y , state &dYdz , const double z);
        void observer(const state& Y, const double z);
        std::vector<double> zs;                                             // Vector container with z values
        std::vector<double> As, dAs, d2A, d3A;                              // Vector container with A in the Einstein frame 
        std::vector<double> lambdas, dlambdas, d2lambda, d3lambda;          // Vector container with lambda = Exp(phi) values
        std::vector<double> Astrings, dAstrings, d2Astring, d3Astring;      // Vector container with A in the string frame
        std::vector<double> Phis, dPhis, d2Phi, d3Phi;                      // Vector container with Dilaton values
        std::vector<double> u0, u2;                                         // Vector container with scalar and tensor glueball potentials
        std::vector<double> aF, bF, cF, l1_2, e2As, e2A;                    // Quantaties useful for spin J fields
        void copy(const IHQCD &ihqcd);
        void finalizeBackground();                                          // Function that slices the std::vector<double> data members
        // Solve the Improved Holographic QCD model
        void solve();
    public:
        IHQCD(const double a0 = 5, const double l0 = 0.0337462, const double z = 10, const double delta = 0.002);
        IHQCD(const IHQCD &ihqcd);
        // Getters of the class
        std::vector<double> getZ();
        std::vector<double> getA();
        std::vector<double> getdA();
        std::vector<double> getd2A();
        std::vector<double> getd3A();
        std::vector<double> getlambda();
        std::vector<double> getdlambda();
        std::vector<double> getd2lambda();
        std::vector<double> getd3lambda();
        std::vector<double> getAstring();
        std::vector<double> getdAstring();
        std::vector<double> getd2Astring();
        std::vector<double> getd3Astring();
        std::vector<double> getPhi();
        std::vector<double> getdPhi();
        std::vector<double> getd2Phi();
        std::vector<double> getd3Phi();
        std::vector<double> getU0();
        std::vector<double> getU2();
        std::vector<double> getaF();
        std::vector<double> getbF();
        std::vector<double> getcF();
        std::vector<double> getl1_2();
        std::vector<double> getE2As();
        std::vector<double> getE2A();
        IHQCD& operator= (const IHQCD &rhs);
        ~IHQCD();
};

IHQCD& ihqcd();
IHQCD& ihqcdU1NNMode();

#endif