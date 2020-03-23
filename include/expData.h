#ifndef EXPDATA_HH
#define EXPDATA_HH

#include <vector>
#include <iostream>
#include <cmath>
#include <utility>
#include <ctime>
#include <cstring>

class ExperimentalData
{
public:
    static constexpr double yrastSpin = 6.5;
    
    //the odd particle's a.m.
    static constexpr double oddSpin_1 = 6.5;
    static constexpr double oddSpin_2 = 4.5;

    //the raw bandhead - units are keV
    // static constexpr double yrastEnergy = 1739.9;
    //transform in MeV
    static constexpr double yrastEnergy = static_cast<double>(1739.9 / 1000.0);
    static constexpr std::size_t nStates_1 = 21;
    static constexpr std::size_t nStates_2 = 17;
    static constexpr std::size_t nStates_3 = 14;
    static constexpr std::size_t nStates_4 = 10;

    //declaring the default constructor and default destructor
public:
    ExperimentalData();
    ~ExperimentalData();

    //methods to initialize the experimental data on class instantation.
public:
    static std::vector<double> init_band1();
    static std::vector<double> init_band2();
    static std::vector<double> init_band3();
    static std::vector<double> init_band4();
    static std::vector<double> generateSpins(double firstSpin, std::size_t length);
    static void arrayPrinter(std::vector<double> &array);
    static void pairPrinter(std::string &bandName, std::vector<double> &spins, std::vector<double> &energies);
    static void generateData(std::vector<double> &spins, std::vector<double> &energies);

    //containers to store the experimental data
public:
    static std::vector<double> spin1;
    static std::vector<double> spin2;
    static std::vector<double> spin3;
    static std::vector<double> spin4;
    static std::vector<double> tsd1;
    static std::vector<double> tsd2;
    static std::vector<double> tsd3;
    static std::vector<double> tsd4;

    //containers to store the unified experimental data
public:
    std::vector<double> spins;
    std::vector<double> energies;
};

#endif // EXPDATA_HH
