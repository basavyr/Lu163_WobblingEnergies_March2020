#ifndef BANDADJUSTER__H
#define BANDADJUSTER__H

#include <vector>
#include <iostream>
#include "../include/expData.h"
#include "../include/energyFormulas.h"

class BandAdjuster
{
public:
    double gamma;
    double beta;
    double I0;
    double V;

public:
    BandAdjuster(double beta, double gamma, double I0, double V);
    ~BandAdjuster();

    //containers to store the necessary numerical data
public:
    static std::vector<double> ExpOld;
    std::vector<double> TSD4Th;
    std::vector<double> ExpNew;
    std::vector<double> Deltas;

    //methods for generating the data
public:
    static std::vector<double> GenerateDeltas(BandAdjuster &obj);
    static std::vector<double> initExpData();
    static std::vector<double> initThData(BandAdjuster &obj);
    static std::vector<double> GenerateNewData(BandAdjuster &obj);
};

#endif // BANDADJUSTER__H
