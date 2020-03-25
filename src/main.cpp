#include "../include/expData.h"
#include "../include/energyFormulas.h"
#include "../include/expData.h"
#include "../include/fitData.h"

int main()
{
    auto nucleus = std::make_unique<ExperimentalData>();
    auto spins = nucleus->spins;
    auto energies = nucleus->energies;
    std::vector<double> TheoreticalEnergies;
    FitData::paramSet bestParams;

    //compute the duration of a process
    auto duration = [](auto start, auto end) {
        auto result = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        return static_cast<double>(result / 1000.0);
    };

    auto startTime = std::chrono::high_resolution_clock::now();
    // FitData::SearchMinimum_RMS(*nucleus, bestParams, 0.38, 17);
    FitData::SearchRMS_SeparateBands(bestParams, 0.38, 17);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "PARAMETERS:"
              << "\n";
    std::cout << bestParams.I0 << " " << bestParams.V << "\n";
    std::cout << "Separate RMS values: \n";
    std::cout << bestParams.E_RMS1 << " " << bestParams.E_RMS2 << " " << bestParams.E_RMS3 << " " << bestParams.E_RMS4 << "\n";
    std::cout << "s= " << bestParams.I0 * bestParams.V << "\n";
    std::cout << "E_RMS =" << FitData::AverageRMS(bestParams) << "\n";
    std::cout << "Process took... " << duration(startTime, endTime) << " s\n";
}