#include "../include/expData.h"
#include "../include/energyFormulas.h"

int main()
{
    auto nucleus = std::make_unique<ExperimentalData>();
    // auto spins = nucleus->spins;
    // auto energies = nucleus->energies;
    EnergyFormulas::moiTuple x(50,0.38,19);
    std::cout << x.A1 << " " << x.A2 << " " << x.A3 << "\n";
}