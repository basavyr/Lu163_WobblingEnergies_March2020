#include "../include/expData.h"
#include "../include/energyFormulas.h"

int main()
{
    auto nucleus = std::make_unique<ExperimentalData>();
    // auto spins = nucleus->spins;
    // auto energies = nucleus->energies;
    // EnergyFormulas::moiTuple x(50, 0.38, 19);
    // std::cout << x.A1 << " " << x.A2 << " " << x.A3 << "\n";
    // for (double i = 0.5; i < 52.5; i += 2.0)
    // {
    //     auto x = new EnergyFormulas::OmegaTuple(i, 6.5, 0.5, 99, 0.38, 19);
    //     std::cout << x->Omega1 << " " << x->Omega2 << "\n";
    // }
    std::cout << EnergyFormulas::TSD1(14.5, 0.5, 99, 0.38, 19) << "\n";
    std::cout << EnergyFormulas::TSD2(14.5, 0.5, 99, 0.38, 19) << "\n";
    std::cout << EnergyFormulas::TSD3(14.5, 0.5, 99, 0.38, 19) << "\n";
    std::cout << EnergyFormulas::TSD4(14.5, 0.5, 99, 0.38, 19) << "\n";
}