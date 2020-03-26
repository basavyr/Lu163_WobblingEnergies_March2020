#include "../include/expData.h"
#include "../include/energyFormulas.h"
#include "../include/expData.h"
#include "../include/fitData.h"
#include "../include/bandAdjuster.h"

void printArray(std::vector<double> &array)
{
    for (int id = 0; id < array.size(); ++id)
    {
        std::cout << array.at(id) << " ";
    }
    std::cout << "\n";
}

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

    //initialize the constructor with the best RMS parameters obtained from the normal fit
    auto adjuster = std::make_unique<BandAdjuster>(0.38, 17, 51.56, 0.0307);
    //create the theoretical data set from the contructor parameters
    adjuster->TSD4Th = BandAdjuster::initThData(*adjuster);
    printArray(adjuster->TSD4Th);

    //based on the values of original experimental data set and the theoretical data set, generate the deviations
    adjuster->Deltas = BandAdjuster::GenerateDeltas(*adjuster);
    printArray(adjuster->Deltas);

    //construct the new experimental data from the deltas
    //the adjusted energies
    adjuster->ExpNew = BandAdjuster::GenerateNewData(*adjuster);
    printArray(adjuster->ExpNew);
    // for (auto &&n : adjuster->ExpNew)
    // std::cout << n << " ";

    //actual fit
    auto startTime = std::chrono::high_resolution_clock::now();
    FitData::SearchRMS_SeparateBands_AdjustedBand(bestParams, *adjuster);
    // // FitData::SearchMinimum_RMS(*nucleus, bestParams, 0.38, 17);
    // FitData::SearchRMS_SeparateBands(bestParams, 0.38, 17);
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "PARAMETERS:"
              << "\n";
    std::cout << bestParams.I0 << " " << bestParams.V << "\n";
    std::cout << "Separate RMS values: \n";
    std::cout << bestParams.E_RMS1 << " " << bestParams.E_RMS2 << " " << bestParams.E_RMS3 << " " << bestParams.E_RMS4 << "\n";
    // std::cout << "s= " << bestParams.I0 * bestParams.V << "\n";
    // std::cout << "E_RMS =" << FitData::AverageRMS(bestParams) << "\n";
    std::cout << "Process took... " << duration(startTime, endTime) << " s\n";
    auto NewAdjustedTheory = FitData::GenerateTheoreticalBand(4, ExperimentalData::spin4, bestParams.I0, bestParams.V, adjuster->beta, adjuster->gamma);
    printArray(NewAdjustedTheory);
}