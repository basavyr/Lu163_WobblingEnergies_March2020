#include "../include/expData.h"

int main()
{
    auto nucleus = std::make_unique<ExperimentalData>();
    auto spins = nucleus->spins;
    auto energies = nucleus->energies;
    std::string name = "all bands";
    ExperimentalData::pairPrinter(name, spins, energies);
}