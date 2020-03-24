#ifndef FITDATA_HH
#define FITDATA_HH

#include <vector>
#include <iostream>
#include <cmath>
#include "../include/expData.h"
#include "energyFormulas.h"

class FitData
{
public:
    FitData(/* args */);
    ~FitData();
    //parameters for storing the best values and limits
public:
    struct paramLimits
    {
        const double I0_left = 48;
        const double I0_right = 58;
        const double I0_step = 0.01;
        const double V_left = 0.1;
        const double V_right = 10;
        const double V_step = 0.01;
    };
    struct paramSet
    {
        double I0;
        double V;
        double E_RMS;
    };

    //methods for calculating the RMS
public:
    static double
    RMS_Calculation(std::vector<double> &expData, std::vector<double> &thData);
    //method to construct the theoretical wobbling excitation container
    static void GenerateTheoreticalEnergies(double i0, double particlePotential, double beta, double gamma, std::vector<double> &energies);
    static void SearchMinimum_RMS(ExperimentalData &obj, FitData::paramSet &bestParams, double beta, double gamma);
};

#endif // FITDATA_HH
