#ifndef GLUONPDF_H
#define GLUONPDF_H

#include "F2.h"
#include "FL.h"
#include "HQCDP.h"

class GluonPDF
{
    public:
        GluonPDF();
        GluonPDF(const GluonPDF &rhs);
        std::vector<double> computePDF(const std::vector<std::vector<double> > & kinPts,
                                       const std::vector<double> pars = {6.46892, -4.69919, 1.12825, 0.664399, -0.0982592, -0.0509552, 0.0176742, -0.0739424, 0.35877});
        GluonPDF& operator=(const GluonPDF &rhs);
        ~GluonPDF();
    private:
        F2 f2;                              // F2 object to be used to compute the PDF
        FL fl;                              // FL object to be used to compute the PDF
        HQCDP pomeron;                      // HQCDP object that will contain the information to compute F2 and FL;
};

#endif