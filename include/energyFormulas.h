#ifndef ENERGYFORMULAS_HH
#define ENERGYFORMULAS_HH

#include <cmath>
#include <vector>
#include <iostream>

class EnergyFormulas
{
public:
    static constexpr double gamma = 17;
    static constexpr double beta = 0.38;
    static constexpr double PI = 3.14159265358979;

public:
    EnergyFormulas(/* args */);
    ~EnergyFormulas();
    //create the necessary tuples for storing the parameters
public:
    struct OmegaTuple
    {
        double Omega1;
        double Omega2;
        OmegaTuple(double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
        {
            Omega1 = EnergyFormulas::omegaFrequency(1, spin, oddSpin, particlePotential, iZero, beta, gamma);
            Omega2 = EnergyFormulas::omegaFrequency(2, spin, oddSpin, particlePotential, iZero, beta, gamma);
        }
    };
    struct moiTuple
    {
        double A1, A2, A3;
        double inertiaFactor(double Ik)
        {
            if (Ik && !isnan(Ik))
                return static_cast<double>(1.0 / (2.0 * Ik));
            return 6969;
        }
        moiTuple()
        {
            /* */
            A1 = 6969;
            A2 = 6969;
            A3 = 6969;
        }
        moiTuple(double iZero, double beta, double gamma)
        {
            auto I1 = EnergyFormulas::rigidMoi(1, iZero, beta, gamma);
            auto I2 = EnergyFormulas::rigidMoi(2, iZero, beta, gamma);
            auto I3 = EnergyFormulas::rigidMoi(3, iZero, beta, gamma);
            A1 = inertiaFactor(I1);
            A2 = inertiaFactor(I2);
            A3 = inertiaFactor(I3);
        }
    };
    //create the helper methods for building the wobbling frequency functions
public:
    static double rigidMoi(int k, double iZero, double beta, double gamma);
    static double H_min(double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma);
    static double B_term(double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma);
    static double C_term(double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma);
    static double omegaFrequency(int n, double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma);
    static double energyExpression(int n1, int n2, double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma);
};

#endif // ENERGYFORMULAS_HH
