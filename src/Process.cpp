#include "Process.h"
#include "methods/vectorOperators.hpp"

void Process::copy(const Process &rhs)
{
    dataPoints = rhs.dataPoints;
}

Process::Process(): dataPoints({}) {}

Process::Process(const Process &proc): dataPoints(proc.dataPoints){}

std::vector<std::vector<double> > Process::getDataPts()
{
    return dataPoints;
}

void Process::setDataPts(const std::vector<std::vector<double> > &pts)
{
    dataPoints = pts ;
}

std::vector<double>  Process::getNeededTVals()
{
    std::vector<double> pts = {0};
    return pts ;
}

double Process::chi2(const std::vector<kinStruct> &Izs, const std::vector<kinStruct> &IzsBar, const std::vector< std::vector < double > > &points)
{
    const std::vector<double> obs = this->diffObsWeighted(Izs, IzsBar, points);
    const int n = obs.size() ;
    double ans = 0.0 ;
    for(int i = 0; i < n ; i++) chi2 += obs[i] * obs[i] ;
    return ans ;
}

Process& Process::operator= (const Process &proc)
{
    if (this == &proc) return *this;
    copy(proc);
    return *this;
}

Process::~Process(){}