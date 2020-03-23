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

//initialize the spins
std::vector<double> ExperimentalData::spin1 = ExperimentalData::generateSpins(ExperimentalData::yrastSpin + 2, ExperimentalData::nStates_1);

std::vector<double> ExperimentalData::spin2 = ExperimentalData::generateSpins(13.5, ExperimentalData::nStates_2);

std::vector<double> ExperimentalData::spin3 = ExperimentalData::generateSpins(16.5, ExperimentalData::nStates_3);

std::vector<double> ExperimentalData::spin4 = ExperimentalData::generateSpins(23.5, ExperimentalData::nStates_4);

//initialize the energies
std::vector<double> ExperimentalData::tsd1 = ExperimentalData::init_band1();

std::vector<double> ExperimentalData::tsd2 = ExperimentalData::init_band2();

std::vector<double> ExperimentalData::tsd3 = ExperimentalData::init_band3();

std::vector<double> ExperimentalData::tsd4 = ExperimentalData::init_band4();

ExperimentalData::ExperimentalData()
{
    ExperimentalData::generateData(spins, energies);
}

ExperimentalData::~ExperimentalData()
{
    int OK = 1;
    std::string firstBand = "TSD1";
    std::string secondBand = "TSD2";
    std::string thirdBand = "TSD3";
    std::string fourthBand = "TSD4";
    // ExperimentalData::pairPrinter(firstBand, spin1, tsd1);
    // ExperimentalData::pairPrinter(secondBand, spin2, tsd2);
    // ExperimentalData::pairPrinter(thirdBand, spin3, tsd3);
    // ExperimentalData::pairPrinter(fourthBand, spin4, tsd4);
    std::cout
        << "Class object finished with status: " << OK;
    std::cout << std::endl;
}

std::vector<double> ExperimentalData::init_band1()
{
    std::vector<double> tsd1 = {1936.5, 2199.6, 2514.5, 2900.8, 3351.1, 3866.4, 4445.0, 5084.0, 5781.0, 6533.6, 7339.1, 8196.9, 9106.6, 10069.2, 11085.7, 12156.8, 13283.0, 14462.3, 15689, 16958, 18262};
    std::vector<double> energies;
    for (auto id = 0; id < tsd1.size(); ++id)
    {
        //raw energy is in keV
        auto currentEnergy = tsd1.at(id);

        //transform in MeV and normalize to bandHead
        currentEnergy = static_cast<double>(currentEnergy / 1000.0 - ExperimentalData::yrastEnergy);
        energies.emplace_back(currentEnergy);
    }
    return energies;
}

std::vector<double> ExperimentalData::init_band2()
{
    std::vector<double> tsd2 = {3079.3, 3486.6, 3958.3, 4492.6, 5088.3, 5742.9, 6454.2, 7220.4, 8040.3, 8913.2, 9839.7, 10819.9, 11854.6, 12943.5, 14086.5, 15284, 16531};
    std::vector<double> energies;
    for (auto id = 0; id < tsd2.size(); ++id)
    {
        //raw energy is in keV
        auto currentEnergy = tsd2.at(id);

        //transform in MeV and normalize to bandHead
        currentEnergy = static_cast<double>(currentEnergy / 1000.0 - ExperimentalData::yrastEnergy);
        energies.emplace_back(currentEnergy);
    }
    return energies;
}
std::vector<double> ExperimentalData::init_band3()
{
    std::vector<double> tsd3 = {3863.6, 4369.2, 4937.2, 5564.2, 6249.3, 6990.5, 7786.4, 8636.2, 9538.7, 10494.5, 11503.7, 12566.7, 13679.1, 14826};
    std::vector<double> energies;
    for (auto id = 0; id < tsd3.size(); ++id)
    {
        //raw energy is in keV
        auto currentEnergy = tsd3.at(id);

        //transform in MeV and normalize to bandHead
        currentEnergy = static_cast<double>(currentEnergy / 1000.0 - ExperimentalData::yrastEnergy);
        energies.emplace_back(currentEnergy);
    }
    return energies;
}
std::vector<double> ExperimentalData::init_band4()
{
    std::vector<double> tsd4 = {6319.9, 6965.0, 7667.2, 8421.8, 9231.8, 10097.2, 11017.7, 11993.4, 13025.0, 14110};
    std::vector<double> energies;
    for (auto id = 0; id < tsd4.size(); ++id)
    {
        //raw energy is in keV
        auto currentEnergy = tsd4.at(id);

        //transform in MeV and normalize to bandHead
        currentEnergy = static_cast<double>(currentEnergy / 1000.0 - ExperimentalData::yrastEnergy);
        energies.emplace_back(currentEnergy);
    }
    return energies;
}

void ExperimentalData::arrayPrinter(std::vector<double> &array)
{
    for (auto &&n : array)
    {
        std::cout << n << " ";
    }
    std::cout << "\n";
}

void ExperimentalData::pairPrinter(std::string &bandName, std::vector<double> &spins, std::vector<double> &energies)
{
    if (spins.size() == energies.size())
    {
        std::cout << bandName << " = {";
        for (int id = 0; id < spins.size(); ++id)
        {
            if (id == spins.size() - 1)
            {
                std::cout << "{ " << spins.at(id) << " , " << energies.at(id) << "}";
            }
            else
            {
                std::cout << "{ " << spins.at(id) << " , " << energies.at(id) << "}, ";
            }
        }
        std::cout << "};"
                  << "\n";
    }
}

void ExperimentalData::generateData(std::vector<double> &spins, std::vector<double> &energies)
{
    for (auto id = 0; id < ExperimentalData::spin1.size(); ++id)
    {
        energies.emplace_back(tsd1.at(id));
        spins.emplace_back(spin1.at(id));
    }
    for (auto id = 0; id < ExperimentalData::spin2.size(); ++id)
    {
        energies.emplace_back(tsd2.at(id));
        spins.emplace_back(spin2.at(id));
    }
    for (auto id = 0; id < ExperimentalData::spin3.size(); ++id)
    {
        energies.emplace_back(tsd3.at(id));
        spins.emplace_back(spin3.at(id));
    }
    for (auto id = 0; id < ExperimentalData::spin4.size(); ++id)
    {
        energies.emplace_back(tsd4.at(id));
        spins.emplace_back(spin4.at(id));
    }
}