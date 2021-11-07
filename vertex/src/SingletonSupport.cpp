// Singletons
#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellId.hpp"

#include "SingletonSupport.hpp"

void SingletonSupport::SetupSingletons(unsigned seed)
{
    SimulationTime::Instance()->SetStartTime(0.0);

    if (seed == 0u)
        RandomNumberGenerator::Instance()->Reseed(time(nullptr));
    else
        RandomNumberGenerator::Instance()->Reseed(seed);

    // //Unnecessary since previous test's tearDown will have cleared:
    // CellPropertyRegistry::Instance()->Clear();

    CellId::ResetMaxCellId();
}

void SingletonSupport::DestroySingletons()
{
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();

    // Destroys properties which are still held by a shared pointer
    CellPropertyRegistry::Instance()->Clear();
}
