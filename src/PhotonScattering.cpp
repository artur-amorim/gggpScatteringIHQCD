#include "PhotonScattering.h"
#include "methods/vectorOperators.hpp"

// Constructor of kinStruct
kinStruct::kinStruct(const double var, const std::vector<double> &izs): kinVar(var), izns(izs) {}

bool kinStruct::operator<(const kinStruct &rhs) const
{
    return kinVar < rhs.kinVar ;
}

bool kinStruct::operator>(const kinStruct &rhs) const
{
    return kinVar > rhs.kinVar ;
}

PhotonScattering::PhotonScattering(): datapts({}){}

PhotonScattering::PhotonScattering(const PhotonScattering &proc):
    datapts(proc.datapts){}

std::vector<std::vector<double> > PhotonScattering::getDataPts()
{
    return datapts;
}

void PhotonScattering::setDataPts(const std::vector<std::vector<double> > &pts)
{
    datapts = pts ;
}

void PhotonScattering::copy(const PhotonScattering &rhs)
{
    datapts = rhs.datapts;
}

void  PhotonScattering::loadData(std::string file_path){}


std::vector<double> PhotonScattering::expVal() 
{
    std::vector<double> pts(0);
    return pts ;
}


std::vector<double> PhotonScattering::expErr()
{
    std::vector<double> pts(0);
    return pts ;
}


std::vector< std::vector<double> >  PhotonScattering::expKinematics()
{
    std::vector< std::vector<double> > pts(0) ;
    return pts ;
}


std::vector<double>  PhotonScattering::getNeededTVals()
{
    std::vector<double> pts = {0};
    return pts ;
}

double PhotonScattering::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    return 0.0 ;
}

double PhotonScattering::IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs)
{
    return 0.0 ;
}


std::vector<kinStruct> PhotonScattering::getIzs(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec)
{
    std::vector<kinStruct> pts;
    return pts ;
}


std::vector<kinStruct>  PhotonScattering::getIzsBar(const std::vector< std::vector<double> >  &points,
                                                     const std::vector<Spectra> &spec, const std::vector<double>  &gs)
{
    std::vector<kinStruct> pts;
    return pts ;
}


std::vector<double>  PhotonScattering::predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                                  const std::vector< std::vector<double> >  &points, const bool savePredictions)
{
    std::vector<double> pts(0) ;
    return pts ;
}


std::vector<double> PhotonScattering::diffObsWeighted(const std::vector<kinStruct> &Izs, const std::vector<kinStruct> &IzsBar, const std::vector< std::vector < double > > &points)
{
    if( points.size() == 0) std::vector< std::vector< double > > points = this->expKinematics() ;   // If points is NULL provide the experimental ones
    const std::vector<double> Opred = this->predict(Izs, IzsBar, points, false) ;                   // Predictions of the model for the process
    const std::vector<double> Oexp  = this->expVal() ;                                              // Experimental values of the process
    const std::vector<double> Oerr  = this->expErr() ;                                              // Experimental errors of the process
    return (Opred - Oexp) / Oerr ;

}

double PhotonScattering::rss(const std::vector<kinStruct> &Izs, const std::vector<kinStruct> &IzsBar, const std::vector< std::vector < double > > &points)
{
    const std::vector<double> obs = this->diffObsWeighted(Izs, IzsBar, points);
    const int n = obs.size() ;
    double chi2 = 0.0 ;
    for(int i = 0; i < n ; i++) chi2 += obs[i] * obs[i] ;
    return chi2 ;
}

PhotonScattering& PhotonScattering::operator= (const PhotonScattering &proc)
{
    if (this == &proc) return *this;
    datapts = proc.datapts;
    return *this;
}

PhotonScattering::~PhotonScattering(){}