#ifndef FITDATA_HH
#define FITDATA_HH

#include <vector>
#include <iostream>
#include <cmath>
#include "../include/expData.h"
#include "../include/energyFormulas.h"
#include "../include/bandAdjuster.h"

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
        const double I0_right = 55;
        const double I0_step = 0.01;
        const double V_left = 0.0001;
        const double V_right = 3;
        const double V_step = 0.001;
    };
    struct paramSet
    {
        double I0;
        double V;
        double E_RMS;
        double E_RMS1, E_RMS2, E_RMS3, E_RMS4;
    };

    //methods for calculating the RMS
public:
    static double
    RMS_Calculation(std::vector<double> &expData, std::vector<double> &thData);
    //method to construct the theoretical wobbling excitation container
    static void GenerateTheoreticalEnergies(double i0, double particlePotential, double beta, double gamma, std::vector<double> &energies);
    static void SearchMinimum_RMS(ExperimentalData &obj, FitData::paramSet &bestParams, double beta, double gamma);
    //search the parameters by getting separate rms for each band
    static void SearchRMS_SeparateBands(FitData::paramSet &bestParams, double beta, double gamma);
    //generates the theoretical band separately (based on each of the experimental ones)
    static std::vector<double> GenerateTheoreticalBand(int bandIndex, std::vector<double> &inputSpins, double iZero, double particlePotential, double beta, double gamma);
    static double AverageRMS(FitData::paramSet &params);

public:
    static void SearchRMS_SeparateBands_AdjustedBand(FitData::paramSet &bestParams, BandAdjuster &obj);
    static void SearchRMS_AdjustedBand(FitData::paramSet &bestParams, BandAdjuster &obj);
};

#endif // FITDATA_HH
