#include "../include/energyFormulas.h"

EnergyFormulas::EnergyFormulas(/* args */)
{
}

EnergyFormulas::~EnergyFormulas()
{
}

double EnergyFormulas::rigidMoi(int k, double iZero, double beta, double gamma)
{
    //work with gamma in radians
    double gammaRad = static_cast<double>(gamma * EnergyFormulas::PI / 180.0);

    const double c1 = sqrt(static_cast<double>(5.0 / (16.0 * EnergyFormulas::PI) * beta));
    const double c2 = sqrt(static_cast<double>(5.0 / (4.0 * EnergyFormulas::PI))) * beta;

    auto frac = static_cast<double>(iZero / (1.0 + c1));
    auto term_k = 1 - c2 * cos(gammaRad + EnergyFormulas::PI * k * 2.0 / 3.0);

    return static_cast<double>(frac * term_k);
}

double EnergyFormulas::H_min(double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
{
    auto I = spin;
    auto j = oddSpin;
    auto V = particlePotential;
    auto gammaRad = gamma * EnergyFormulas::PI / 180.0;
    EnergyFormulas::moiTuple mois(iZero, beta, gamma);
    const double pi6 = static_cast<double>(EnergyFormulas::PI / 6.0);
    auto trigFunction = sin(gammaRad + pi6);
    auto term1 = mois.A2 + mois.A3 * (I + j) * 0.5;
    auto term2 = mois.A1 * pow(I - j, 2);
    auto term3 = V * (2.0 * j - 1.0) / (j + 1.0) * trigFunction;
    auto result = term1 + term2 - term3;
    if (!isnan(result))
        return result;
    return 6969;
}

double EnergyFormulas::B_term(double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
{
    //declaring the constants and pre-terms
    auto I = spin;
    auto j = oddSpin;
    auto V = particlePotential;
    auto gammaRad = gamma * EnergyFormulas::PI / 180.0;
    EnergyFormulas::moiTuple mois(iZero, beta, gamma);
    auto sing = sin(gammaRad);
    auto cosg = cos(gammaRad);

    //declaring the actual terms of the function
    auto term1 = (2.0 * I - 1.0) * (mois.A3 - mois.A1) + 2.0 * j * mois.A1;
    auto term2 = (2.0 * I - 1.0) * (mois.A2 - mois.A1) + 2.0 * j * mois.A1;
    auto term3 = 8.0 * mois.A2 * mois.A3 * I * j;
    auto term4 = (2.0 * j - 1.0) * (mois.A3 - mois.A1) + 2.0 * I * j * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * sqrt(3.0) * (sqrt(3.0) * cosg + sing);
    auto term5 = (2.0 * j - 1.0) * (mois.A2 - mois.A1) + 2.0 * I * j * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * 2.0 * sqrt(3.0) * sing;
    auto result = term1 * term2 + term3 + term4 * term5;
    if (!isnan(result))
        return result;
    return 6969;
}

double EnergyFormulas::C_term(double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
{
    //declaring the constants and pre-terms
    auto I = spin;
    auto j = oddSpin;
    auto V = particlePotential;
    auto gammaRad = gamma * EnergyFormulas::PI / 180.0;
    EnergyFormulas::moiTuple mois(iZero, beta, gamma);
    auto sing = sin(gammaRad);
    auto cosg = cos(gammaRad);

    //declaring the actual terms of the function
    auto term1 = (2.0 * I - 1.0) * (mois.A3 - mois.A1) + 2.0 * j * mois.A1;
    auto term2 = (2.0 * j - 1.0) * (mois.A3 - mois.A1) + 2.0 * I * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * sqrt(3.0) * (sqrt(3.0) * cosg + sing);
    auto bracket1 = term1 * term2 - 4.0 * I * j * pow(mois.A3, 2);
    auto term3 = (2.0 * I - 1.0) * (mois.A2 - mois.A1) + 2.0 * j * mois.A1;
    auto term4 = (2.0 * j - 1.0) * (mois.A3 - mois.A1) + 2.0 * I * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * 2.0 * sqrt(3.0) * sing;
    auto bracket2 = term3 * term4 - 4.0 * I * j * pow(mois.A2, 2);
    auto result = bracket1 * bracket2;
    if (!isnan(result))
        return result;
    return 6969;
}

double EnergyFormulas::omegaFrequency(int n, double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
{
    auto B = EnergyFormulas::B_term(spin, oddSpin, particlePotential, iZero, beta, gamma);
    if (isnan(B))
        return 6969;
    auto C = EnergyFormulas::C_term(spin, oddSpin, particlePotential, iZero, beta, gamma);
    if (isnan(C))
        return 6969;
    auto root = sqrt(pow(B, 2) - 4.0 * C);
    if (isnan(root))
        return 6969;
    double result = 0;
    if (n == 1)
    {
        result = 0.5 * (-B - root);
    }
    else
    {
        result = 0.5 * (-B + root);
    }
    if (!isnan(result))
        return sqrt(result);
    return 6969;
}
double EnergyFormulas::energyExpression(int n1, int n2, double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
{
    auto I = spin;
    auto j = oddSpin;
    auto V = particlePotential;
    auto minEnergy = EnergyFormulas::H_min(I, j, V, iZero, beta, gamma);
    //continue only if the minimum term is real
    if (isnan(minEnergy))
        return 6969;
    OmegaTuple Omegas(I, j, V, iZero, beta, gamma);
    if (Omegas.Omega1 == 6969 || Omegas.Omega2 == 6969)
        return 6969;
    auto Phonon_One = (n1 + 0.5) * Omegas.Omega1;
    auto Phonon_Two = (n2 + .05) * Omegas.Omega2;
    auto energy_I = minEnergy + Phonon_One + Phonon_Two;
    if (!isnan(energy_I))
        return energy_I;
    return 6969;
}