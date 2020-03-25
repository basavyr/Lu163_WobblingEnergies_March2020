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
    auto result = static_cast<double>(sqrt(sum / (size - 1.0)));
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

std::vector<double> FitData::GenerateTheoreticalBand(int bandIndex, std::vector<double> &inputSpins, double iZero, double particlePotential, double beta, double gamma)
{
    std::vector<double> results;
    auto V = particlePotential;
    auto I0 = iZero;
    double currentEnergy = -1.0;
    for (auto id = 0; id < inputSpins.size(); ++id)
    {
        auto I = inputSpins.at(id);
        switch (bandIndex)
        {
        case 1:
            currentEnergy = EnergyFormulas::TSD1(I, V, I0, beta, gamma);
            break;
        case 2:
            currentEnergy = EnergyFormulas::TSD2(I, V, I0, beta, gamma);
            break;
        case 3:
            currentEnergy = EnergyFormulas::TSD3(I, V, I0, beta, gamma);
            break;
        case 4:
            currentEnergy = EnergyFormulas::TSD4(I, V, I0, beta, gamma);
            break;
        }
        if (currentEnergy != 6969 && currentEnergy >= 0)
            results.emplace_back(currentEnergy);
    }
    return results;
}

void FitData::SearchRMS_SeparateBands(FitData::paramSet &bestParams, double beta, double gamma)
{
    double minValue_band1 = 987654321.0;
    double minValue_band2 = 987654321.0;
    double minValue_band3 = 987654321.0;
    double minValue_band4 = 987654321.0;
    FitData::paramLimits limits;
    for (double I0 = limits.I0_left; I0 <= limits.I0_right; I0 += limits.I0_step)
    {
        for (double V = limits.V_left; V <= limits.V_right; V += limits.V_step)
        {
            std::vector<double> band1 = FitData::GenerateTheoreticalBand(1, ExperimentalData::spin1, I0, V, beta, gamma);
            std::vector<double> band2 = FitData::GenerateTheoreticalBand(2, ExperimentalData::spin2, I0, V, beta, gamma);
            std::vector<double> band3 = FitData::GenerateTheoreticalBand(3, ExperimentalData::spin3, I0, V, beta, gamma);
            std::vector<double> band4 = FitData::GenerateTheoreticalBand(4, ExperimentalData::spin4, I0, V, beta, gamma);
            auto RMS1 = FitData::RMS_Calculation(ExperimentalData::tsd1, band1);
            auto RMS2 = FitData::RMS_Calculation(ExperimentalData::tsd2, band2);
            auto RMS3 = FitData::RMS_Calculation(ExperimentalData::tsd3, band3);
            auto RMS4 = FitData::RMS_Calculation(ExperimentalData::tsd4, band4);
            if ((RMS1 != 6969 && RMS1 <= minValue_band1) && (RMS2 != 6969 && RMS2 <= minValue_band2) && (RMS3 != 6969 && RMS3 <= minValue_band3) && (RMS4 != 6969 && RMS4 <= minValue_band4))
            {
                minValue_band1 = RMS1;
                minValue_band2 = RMS2;
                minValue_band3 = RMS3;
                minValue_band4 = RMS4;
                bestParams.I0 = I0;
                bestParams.V = V;
                bestParams.E_RMS1 = minValue_band1;
                bestParams.E_RMS2 = minValue_band2;
                bestParams.E_RMS3 = minValue_band3;
                bestParams.E_RMS4 = minValue_band4;
            }
        }
    }
}

double FitData::AverageRMS(FitData::paramSet &params)
{
    //only work if the params are real
    auto rms1 = params.E_RMS1;
    auto rms2 = params.E_RMS2;
    auto rms3 = params.E_RMS3;
    auto rms4 = params.E_RMS4;
    double avg = 6969;
    if (rms1 != 6969 && rms2 != 6969 && rms3 != 6969 && rms4 != 6969)
    {
        avg = 0.25 * (rms1 + rms2 + rms3 + rms4);
    }
    if (avg != 6969)
        return avg;
    return 0;
}