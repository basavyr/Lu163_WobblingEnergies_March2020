#ifndef FITDATA_HH
#define FITDATA_HH

#include <vector>
#include <iostream>

class FitData
{
public:
    FitData(/* args */);
    ~FitData();

    //methods for calculating the RMS
public:
    static double RMS_Calculation(std::vector<double> &expData, std::vector<double> &thData);
    //method to construct the theoretical wobbling excitation container
    static void GenerateTheoreticalEnergies(std::vector<double> &energies);
};

#endif // FITDATA_HH
