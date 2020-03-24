#include "../include/fitData.h"
FitData::FitData(/* args */)
{
}

FitData::~FitData()
{
}

double FitData::RMS_Calculation(std::vector<double> &expData, std::vector<double> &thData)
{
    //stop if the arrays are not equal in size
    if (expData.size() != thData.size())
    {
        // std::cout << "AND IT ACTUALLY FAILED..." << "\n";
        return 6969;
    }
    double sum = 0;
    std::size_t size = expData.size();
    int errorCheck = 0;
    int ok = 1;
    for (auto id = 0; id < size && ok; ++id)
    {
        auto currentRMS = pow(expData.at(id) - thData.at(id), 2);
        if (!isnan(currentRMS))
        {
            sum += currentRMS;
            errorCheck++;
        }
        else
        {
            ok = 0;
            id = size + 1;
        }
    }
    auto result = static_cast<double>(sqrt(sum / (size + 1.0)));
    if (!isnan(result) && errorCheck == size && ok)
    {
        return result;
    }
    return 6969;
}

void FitData::GenerateTheoreticalEnergies(double i0, double particlePotential, double beta, double gamma, std::vector<double> &energies)
{
    // //continue only if the array is empty
    // if (energies.size())
    //     return;
    auto V = particlePotential;
    int errorChecker = 0;
    //tsd1
    for (auto id = 0; id < ExperimentalData::spin1.size(); ++id)
    {
        auto I = ExperimentalData::spin1.at(id);
        auto currentEnergy = EnergyFormulas::TSD1(I, V, i0, beta, gamma);
        if (currentEnergy != 6969)
        {
            energies.emplace_back(currentEnergy);
            errorChecker++;
        }
    }
    //tsd2
    for (auto id = 0; id < ExperimentalData::spin2.size(); ++id)
    {
        auto I = ExperimentalData::spin2.at(id);
        auto currentEnergy = EnergyFormulas::TSD2(I, V, i0, beta, gamma);
        if (currentEnergy != 6969)
        {
            energies.emplace_back(currentEnergy);
            errorChecker++;
        }
    }
    //tsd3
    for (auto id = 0; id < ExperimentalData::spin3.size(); ++id)
    {
        auto I = ExperimentalData::spin3.at(id);
        auto currentEnergy = EnergyFormulas::TSD3(I, V, i0, beta, gamma);
        if (currentEnergy != 6969)
        {
            energies.emplace_back(currentEnergy);
            errorChecker++;
        }
    }
    //tsd4
    for (auto id = 0; id < ExperimentalData::spin4.size(); ++id)
    {
        auto I = ExperimentalData::spin4.at(id);
        auto currentEnergy = EnergyFormulas::TSD4(I, V, i0, beta, gamma);
        if (currentEnergy != 6969)
        {
            energies.emplace_back(currentEnergy);
            errorChecker++;
        }
    }
    // if (errorChecker != (ExperimentalData::spin1.size() + ExperimentalData::spin2.size() + ExperimentalData::spin3.size() + ExperimentalData::spin4.size()))
    // std::cout << "SHOULD FAIL IN MIN_RMS CALCULATION FOR I0= " << i0 <<" AND V= " << V  << "\n";
}

void FitData::SearchMinimum_RMS(ExperimentalData &obj, FitData::paramSet &bestParams, double beta, double gamma)
{
    double minValue = 98765432101.23456789;
    auto ExperimentalEnergies = obj.energies;
    FitData::paramLimits limits;
    //start the searching procedure
    for (double I0 = limits.I0_left; I0 <= limits.I0_right; I0 += limits.I0_step)
    {
        for (double V = limits.V_left; V <= limits.V_right; V += limits.V_step)
        {
            std::vector<double> TheoreticalEnergies;
            FitData::GenerateTheoreticalEnergies(I0, V, beta, gamma, TheoreticalEnergies);
            auto currentRMS = FitData::RMS_Calculation(ExperimentalEnergies, TheoreticalEnergies);
            if (currentRMS != 6969 && currentRMS <= minValue)
            {
                minValue = currentRMS;
                bestParams.I0 = I0;
                bestParams.V = V;
                bestParams.E_RMS = minValue;
            }
        }
    }
}