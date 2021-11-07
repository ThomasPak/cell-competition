#include "AbstractPhaseBasedCellCycleModel.hpp"

#include "DeathClockCellKiller.hpp"

template<unsigned DIM>
DeathClockCellKiller<DIM>::DeathClockCellKiller(
        AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template<unsigned DIM>
void DeathClockCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (auto cell_iter = this->mpCellPopulation->Begin();
        cell_iter != this->mpCellPopulation->End();
        ++cell_iter)
    {
        // If cell is apoptotic, skip cell
        if (cell_iter->HasApoptosisBegun())
            continue;

        // If cell cycle model is phase-based, kill only cells in G1 phase
        if (AbstractPhaseBasedCellCycleModel *p_ccm =
                dynamic_cast<AbstractPhaseBasedCellCycleModel*>(cell_iter->GetCellCycleModel()))
        {
            // If cell is not in G1 phase, skip cell
            if (p_ccm->GetCurrentCellCyclePhase() != G_ONE_PHASE)
                continue;
        }

        double tau = cell_iter->GetCellData()->GetItem("tau");
        double death_threshold =
            cell_iter->GetCellData()->GetItem("death threshold");

        // Sanity checks
        assert(tau >= 0.0);
        assert(death_threshold > 0.0);

        // Start apoptosis if tau has passed death threshold
        if (tau >= death_threshold)
        {
            cell_iter->StartApoptosis(false);
        }
    }
}

template<unsigned DIM>
void DeathClockCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}
// Explicit instantiation
template class DeathClockCellKiller<1>;
template class DeathClockCellKiller<2>;
template class DeathClockCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeathClockCellKiller)
