#ifndef PROCESS_H
#define PROCESS_H

#include <vector>
#include "Spectra.h"

// Auxiliary struct to reduce the number of IzN and IzNBarcomputations
struct kinStruct
{
    double kinVar;
    std::vector<double> izns;
    kinStruct(const double var = 0, const std::vector<double> &izs = {});
    bool operator<(const kinStruct &rhs) const;
    bool operator>(const kinStruct &rhs) const;
};

class Process
{
    private:
        std::vector<std::vector<double> > dataPoints;
    protected:
        virtual void copy(const Process &proc);
    public:
        Process();
        Process(const Process &proc);
        std::vector< std::vector<double> > getDataPts();
        void setDataPts(const std::vector< std::vector<double> > &pts );
        virtual std::vector<double> expVal() = 0;
        virtual std::vector<double>  expErr() = 0;
        virtual std::vector<std::vector<double> >  expKinematics() = 0;
        std::vector<double>  getNeededTVals();
        virtual double IzN(const std::vector<double> &kin, const Reggeon &reg) = 0;
        virtual double IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs) = 0;
        virtual std::vector<kinStruct>  getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec) = 0;
        virtual std::vector<kinStruct>  getIzsBar(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec,
                                                  const std::vector<double>  &gs) = 0;
        virtual std::vector<double>  predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                            const std::vector< std::vector<double> >  &points, const bool savePredictions = false) = 0;
        virtual std::vector<double>  diffObsWeighted(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar, 
                                                    const std::vector< std::vector<double> >  &points) = 0;
        double chi2(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar, const std::vector< std::vector<double> >  &points) ;   
        Process& operator= (const Process &proc);
        virtual ~Process();
};

#endif