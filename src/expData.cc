#include "../include/expData.h"

std::vector<double> ExperimentalData::generateSpins(double firstSpin, std::size_t length)
{
    std::vector<double> spins;
    for (auto id = 0; id < length; ++id)
        spins.emplace_back(static_cast<double>(firstSpin + id * 2.0));
    if (spins.size())
        return spins;
    return spins;
}

std::vector<double> ExperimentalData::spin1 = ExperimentalData::generateSpins(ExperimentalData::yrastSpin, ExperimentalData::nStates_1);

std::vector<double> ExperimentalData::spin2 = ExperimentalData::generateSpins(13.5, ExperimentalData::nStates_2);

std::vector<double> ExperimentalData::spin3 = ExperimentalData::generateSpins(16.5, ExperimentalData::nStates_3);

std::vector<double> ExperimentalData::spin4 = ExperimentalData::generateSpins(23.5, ExperimentalData::nStates_4);

ExperimentalData::ExperimentalData()
{
}

ExperimentalData::~ExperimentalData()
{
    int OK = 1;
    std::cout << "Class object finished with status: " << OK;
    std::cout << std::endl;
}

// std::vector<double> ExperimentalData::init_band1()
// {
//     std::vector<double> spins;
//     std::vector<double> energies;
//     for (auto id = 0; id < spins.size(); ++id)
//     {
//     }
// }

// std::vector<double> ExperimentalData::init_band2()
// {
// }
// std::vector<double> ExperimentalData::init_band3()
// {
// }
// std::vector<double> ExperimentalData::init_band4()
// {
// }
