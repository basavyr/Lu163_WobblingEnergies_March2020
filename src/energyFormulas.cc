#include "../include/energyFormulas.h"
#include "../include/expData.h"

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

    const double c1 = sqrt(static_cast<double>(5.0 / (16.0 * EnergyFormulas::PI)));
    const double c2 = sqrt(static_cast<double>(5.0 / (4.0 * EnergyFormulas::PI))) * beta;

    auto frac = static_cast<double>(iZero / (1.0 + c1 * beta));
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

    //stop if the MOIS are not physical
    if (mois.A1 == 6969 || mois.A2 == 6969 || mois.A3 == 6969)
        return 6969;

    const double pi6 = static_cast<double>(EnergyFormulas::PI / 6.0);
    auto trigFunction = sin(gammaRad + pi6);
    auto term1 = (mois.A2 + mois.A3) * (I + j) * 0.5;
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

    //stop if the MOIS are not physical
    if (mois.A1 == 6969 || mois.A2 == 6969 || mois.A3 == 6969)
        return 6969;

    auto sing = sin(gammaRad);
    auto cosg = cos(gammaRad);

    //declaring the actual terms of the function
    auto term1 = (2.0 * I - 1.0) * (mois.A3 - mois.A1) + 2.0 * j * mois.A1;
    auto term2 = (2.0 * I - 1.0) * (mois.A2 - mois.A1) + 2.0 * j * mois.A1;
    auto term3 = 8.0 * mois.A2 * mois.A3 * I * j;
    auto term4 = (2.0 * j - 1.0) * (mois.A3 - mois.A1) + 2.0 * I * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * sqrt(3.0) * (sqrt(3.0) * cosg + sing);
    auto term5 = (2.0 * j - 1.0) * (mois.A2 - mois.A1) + 2.0 * I * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * 2.0 * sqrt(3.0) * sing;
    auto result = -1.0 * (term1 * term2 + term3 + term4 * term5);
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

    //stop if the MOIS are not physical
    if (mois.A1 == 6969 || mois.A2 == 6969 || mois.A3 == 6969)
        return 6969;

    //declaring the actual terms of the function
    auto term1 = (2.0 * I - 1.0) * (mois.A3 - mois.A1) + 2.0 * j * mois.A1;
    auto term2 = (2.0 * j - 1.0) * (mois.A3 - mois.A1) + 2.0 * I * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * sqrt(3.0) * (sqrt(3.0) * cosg + sing);
    auto bracket1 = term1 * term2 - 4.0 * I * j * pow(mois.A3, 2);
    auto term3 = (2.0 * I - 1.0) * (mois.A2 - mois.A1) + 2.0 * j * mois.A1;
    auto term4 = (2.0 * j - 1.0) * (mois.A2 - mois.A1) + 2.0 * I * mois.A1 + V * (2.0 * j - 1.0) / (j * (j + 1.0)) * 2.0 * sqrt(3.0) * sing;
    auto bracket2 = term3 * term4 - 4.0 * I * j * pow(mois.A2, 2);
    auto result = bracket1 * bracket2;
    if (!isnan(result))
        return result;
    return 6969;
}

double EnergyFormulas::omegaFrequency(int n, double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
{
    auto B = EnergyFormulas::B_term(spin, oddSpin, particlePotential, iZero, beta, gamma);

    //***************
    //debug
    // std::cout << "B= " << B << "\n";
    //***************

    if (isnan(B))
        return 6969;

    auto C = EnergyFormulas::C_term(spin, oddSpin, particlePotential, iZero, beta, gamma);

    //***************
    //debug
    // std::cout << "C= " << C << "\n";
    //***************

    if (isnan(C))
        return 6969;

    auto root = sqrt(pow(B, 2) - 4.0 * C);
    //***************
    //debug
    // std::cout << "root= " << root << "\n";
    //***************

    if (isnan(root))
        return 6969;

    double result = -1.0;
    if (n == 1)
    {
        result = 0.5 * (-B - root);
        //***************
        //debug
        // std::cout << "Omega1= " << result << "\n";
        //***************
    }
    else
    {
        result = 0.5 * (-B + root);
        //***************
        //debug
        // std::cout << "Omega2= " << result << "\n";
        //***************
    }
    if (!isnan(result) && result > 0.0)
    {
        // std::cout << "OMEGAS ARE OK " << sqrt(result) << "\n";
        return sqrt(result);
    }
    return 6969;
}
double EnergyFormulas::energyExpression(int n1, int n2, double spin, double oddSpin, double particlePotential, double iZero, double beta, double gamma)
{
    auto I = spin;
    //***************
    //debug
    // std::cout << "I= " << I << "\n";
    //***************
    auto j = oddSpin;
    //***************
    //debug
    // std::cout << "j= " << j << "\n";
    //***************
    auto V = particlePotential;
    //***************
    //debug
    // std::cout << "V= " << V << "\n";
    //***************
    auto minEnergy = EnergyFormulas::H_min(I, j, V, iZero, beta, gamma);
    //***************
    //debug
    // std::cout << "H_min= " << minEnergy << "\n";
    //***************

    //continue only if the minimum term is real
    if (isnan(minEnergy))
        return 6969;

    OmegaTuple Omegas(I, j, V, iZero, beta, gamma);
    //***************
    //debug
    // std::cout << "Omegas= " << Omegas.Omega1 << " " << Omegas.Omega2 << "\n";
    //***************
    if (Omegas.Omega1 == 6969 || Omegas.Omega2 == 6969)
        return 6969;

    auto Phonon_One = (n1 + 0.5) * Omegas.Omega1;
    //***************
    //debug
    // std::cout << "PhononOne= " << Phonon_One << "\n";
    //***************
    auto Phonon_Two = (n2 + 0.5) * Omegas.Omega2;
    //***************
    //debug
    // std::cout << "PhononOne= " << Phonon_Two << "\n";
    //***************

    auto energy_I = minEnergy + Phonon_One + Phonon_Two;
    //***************
    //debug
    // std::cout << "energy_I= " << energy_I << "\n";
    //***************
    if (!isnan(energy_I))
        return energy_I;
    return 6969;
}

double EnergyFormulas::TSD1(double spin, double particlePotential, double iZero, double beta, double gamma)
{
    auto V = particlePotential;
    //***************
    //debug
    // std::cout << "V= " << V << "\n";
    //***************
    auto I = spin;
    //***************
    //debug
    // std::cout << "I= " << I << "\n";
    //***************
    auto I_0 = iZero;
    //***************
    //debug
    // std::cout << "I_0= " << I_0 << "\n";
    //***************
    auto j = ExperimentalData::oddSpin_1;
    //***************
    //debug
    // std::cout << "j= " << j << "\n";
    //***************

    //THE BANDHEAD ENERGY (state with I=13/2)
    double yrastEnergy = energyExpression(0, 0, ExperimentalData::yrastSpin, ExperimentalData::oddSpin_1, V, I_0, beta, gamma);
    //***************
    //debug
    // std::cout << "yrastEnergy= " << yrastEnergy << "\n";
    //***************

    //stop if the yrast energy is not physical
    if (yrastEnergy == 6969 || isnan(yrastEnergy))
        return 6969;

    auto energy_I = EnergyFormulas::energyExpression(0, 0, I, j, V, I_0, beta, gamma);
    //***************
    //debug
    // std::cout << "energy_I= " << energy_I << "\n";
    //***************

    //stop of the current wobbling energy is not real
    if (energy_I == 6969)
        return 6969;

    //calculate the excitation energy
    auto result = energy_I - yrastEnergy;
    //***************
    //debug
    // std::cout << "E_exc= " << result << "\n";
    //***************

    if (!isnan(result) && result >= 0.0)
        return result;
    return 6969;
}

double EnergyFormulas::TSD2(double spin, double particlePotential, double iZero, double beta, double gamma)
{
    auto V = particlePotential;
    auto I = spin - 1.0;
    auto I_0 = iZero;
    auto j = ExperimentalData::oddSpin_1;

    //THE BANDHEAD ENERGY (state with I=13/2)
    double yrastEnergy = energyExpression(0, 0, ExperimentalData::yrastSpin, ExperimentalData::oddSpin_1, V, I_0, beta, gamma);

    //stop if the yrast energy is not physical
    if (yrastEnergy == 6969 || isnan(yrastEnergy))
        return 6969;

    auto energy_I = EnergyFormulas::energyExpression(1, 0, I, j, V, I_0, beta, gamma);

    //stop of the current wobbling energy is not real
    if (energy_I == 6969)
        return 6969;

    //calculate the excitation energy
    auto result = energy_I - yrastEnergy;
    if (!isnan(result) && result >= 0.0)
        return result;
    return 6969;
}

double EnergyFormulas::TSD3(double spin, double particlePotential, double iZero, double beta, double gamma)
{
    auto V = particlePotential;
    auto I = spin - 2.0;
    auto I_0 = iZero;
    auto j = ExperimentalData::oddSpin_1;

    //THE BANDHEAD ENERGY (state with I=13/2)
    double yrastEnergy = energyExpression(0, 0, ExperimentalData::yrastSpin, ExperimentalData::oddSpin_1, V, I_0, beta, gamma);

    //stop if the yrast energy is not physical
    if (yrastEnergy == 6969 || isnan(yrastEnergy))
        return 6969;

    auto energy_I = EnergyFormulas::energyExpression(2, 0, I, j, V, I_0, beta, gamma);

    //stop of the current wobbling energy is not real
    if (energy_I == 6969)
        return 6969;

    //calculate the excitation energy
    auto result = energy_I - yrastEnergy;
    if (!isnan(result) && result >= 0.0)
        return result;
    return 6969;
}

double EnergyFormulas::TSD4(double spin, double particlePotential, double iZero, double beta, double gamma)
{
    auto V = particlePotential;
    auto I = spin - 3.0;
    auto I_0 = iZero;
    auto j = ExperimentalData::oddSpin_2;

    //THE BANDHEAD ENERGY (state with I=13/2)
    double yrastEnergy = energyExpression(0, 0, ExperimentalData::yrastSpin, ExperimentalData::oddSpin_1, V, I_0, beta, gamma);

    //stop if the yrast energy is not physical
    if (yrastEnergy == 6969 || isnan(yrastEnergy))
        return 6969;

    auto energy_I = EnergyFormulas::energyExpression(3, 0, I, j, V, I_0, beta, gamma);

    //stop of the current wobbling energy is not real
    if (energy_I == 6969)
        return 6969;

    //calculate the excitation energy
    auto result = energy_I - yrastEnergy;
    if (!isnan(result) && result >= 0.0)
        return result;
    return 6969;
}