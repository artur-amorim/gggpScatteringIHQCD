#ifndef GAMMASCAT_H
#define GAMMASCAT_H

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

class PhotonScattering
{
    private :
        std::vector< std::vector<double> > datapts ;
    protected:
        virtual void copy(const PhotonScattering & rhs);
        virtual void  loadData(std::string file_path = "");
    public :
        PhotonScattering();
        PhotonScattering(const PhotonScattering &proc);
        std::vector< std::vector<double> > getDataPts();
        void setDataPts(const std::vector< std::vector<double> > &pts );
        virtual std::vector<double> expVal();
        virtual std::vector<double>  expErr();
        virtual std::vector<std::vector<double> >  expKinematics();
        std::vector<double>  getNeededTVals();
        virtual double IzN(const std::vector<double> &kin, const Reggeon &reg);
        virtual double IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs);
        virtual std::vector<kinStruct>  getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec);
        virtual std::vector<kinStruct>  getIzsBar(const std::vector< std::vector<double> >  &points,
                                                  const std::vector<Spectra> &spec,
                                                  const std::vector<double>  &gs);
        virtual std::vector<double>  predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                                const std::vector< std::vector<double> >  &points, const bool savePredictions = false);
        virtual std::vector<double>  diffObsWeighted(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar, const std::vector< std::vector<double> >  &points);
        double rss(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar, const std::vector< std::vector<double> >  &points) ;   
        PhotonScattering& operator= (const PhotonScattering &procs);
        virtual ~PhotonScattering();
};
#endif