#include "../include/bandAdjuster.h"

BandAdjuster::BandAdjuster(double b1, double g1, double I01, double V1)
{
    gamma = g1;
    beta = b1;
    I0 = I01;
    V = V1;
}

BandAdjuster::~BandAdjuster()
{
}

std::vector<double> BandAdjuster::ExpOld = BandAdjuster::initExpData();

std::vector<double> BandAdjuster::initExpData()
{
    //the experimental values for the energies in the tsd4 band
    std::vector<double> result;
    for (int id = 0; id < ExperimentalData::tsd4.size(); ++id)
    {
        auto currentEnergy = ExperimentalData::tsd4.at(id);
        result.emplace_back(currentEnergy);
    }
    return result;
}

std::vector<double> BandAdjuster::initThData(BandAdjuster &obj)
{
    //the theoretical values for the energies in the tsd4 band
    //obtain with the classic algorithm -> best possible rms this is
    std::vector<double> result;
    auto beta = obj.beta;
    auto gamma = obj.gamma;
    auto V = obj.V;
    auto I0 = obj.I0;
    for (int id = 0; id < ExperimentalData::tsd4.size(); ++id)
    {
        auto I = ExperimentalData::spin4.at(id);
        auto currentEnergy = EnergyFormulas::TSD4(I, V, I0, beta, gamma);
        result.emplace_back(currentEnergy);
    }
    return result;
}

std::vector<double> BandAdjuster::GenerateDeltas(BandAdjuster &obj)
{
    std::vector<double> deltas;
    for (int id = 0; id < obj.TSD4Th.size(); ++id)
    {
        auto expval = BandAdjuster::ExpOld.at(id);
        auto thval = obj.TSD4Th.at(id);
        deltas.emplace_back(thval - expval);
    }
    return deltas;
}

std::vector<double> BandAdjuster::GenerateNewData(BandAdjuster &obj)
{
    std::vector<double> result;
    for (int id = 0; id < obj.ExpOld.size(); ++id)
    {
        auto currentValue = BandAdjuster::ExpOld.at(id) - obj.Deltas.at(id);
        result.emplace_back(currentValue);
    }
    return result;
}