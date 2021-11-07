#include "DeathClockTargetAreaModifier.hpp"
#include "ApoptoticCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"

template<unsigned DIM>
DeathClockTargetAreaModifier<DIM>::DeathClockTargetAreaModifier()
    : AbstractTargetAreaModifier<DIM>(),
      mAgeToStartGrowing(DOUBLE_UNSET),
      mGrowthRate(DOUBLE_UNSET)
{
}

template<unsigned DIM>
DeathClockTargetAreaModifier<DIM>::~DeathClockTargetAreaModifier()
{
}

template<unsigned DIM>
void DeathClockTargetAreaModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area =
        pCell->GetCellData()->GetItem("reference target area");

    // The target area of an apoptotic cell decreases linearly to zero
    if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        ///\todo: which cells are apoptotic? if they get apoptotic during G2-phase then this line has to be changed
        double time_spent_apoptotic = SimulationTime::Instance()->GetTime() - pCell->GetStartOfApoptosisTime();
        cell_target_area = cell_target_area - 0.5*cell_target_area/(pCell->GetApoptosisTime())*time_spent_apoptotic;

        // Don't allow a negative target area
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else
    {
        bool cell_is_differentiated = pCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>();
        if (!cell_is_differentiated)
        {
            /*
             * If not using a phase-based cell-cycle model, the target area of a proliferating cell increases
             * linearly from mReferenceTargetArea as soon as the cell's age exceeds mAgeToStartGrowing, with
             * growth rate mGrowthRate.
             */
            double growth_start_time = mAgeToStartGrowing;
            double growth_rate = mGrowthRate;
            if (growth_start_time == DOUBLE_UNSET)
            {
                if (dynamic_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel()) == nullptr)
                {
                    EXCEPTION("If SetAgeToStartGrowing() has not been called, a subclass of AbstractPhaseBasedCellCycleModel must be used");
                }

                /*
                 * If using a phase-based cell-cycle model, the target area of a proliferating cell increases
                 * linearly from mReferenceTargetArea to 2*mReferenceTargetArea during the cell's G2 phase.
                 */
                AbstractPhaseBasedCellCycleModel* p_model = static_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel());
                growth_start_time = p_model->GetMDuration() + p_model->GetG1Duration() + p_model->GetSDuration();

                double g2_duration = p_model->GetG2Duration();
                growth_rate = cell_target_area/g2_duration;
            }
            else if (mGrowthRate == DOUBLE_UNSET)
            {
                EXCEPTION("If SetAgeToStartGrowing() has been called, then SetGrowthRate() must also be called");
            }

            double time_spent_growing = pCell->GetAge() - growth_start_time;

            // The target area of a proliferating cell increases linearly
            if (time_spent_growing > 0)
            {
                cell_target_area += time_spent_growing*growth_rate;
            }

            /**
             * At division, daughter cells inherit the cell data array from the mother cell.
             * Here, we assign the target area that we want daughter cells to have to cells
             * that we know to divide in this time step.
             *
             * \todo This is a little hack that we might want to clean up in the future.
             */
            if (pCell->ReadyToDivide())
            {
                cell_target_area = this->mReferenceTargetArea;
            }
        }
    }

    // Set cell data
    pCell->GetCellData()->SetItem("target area", cell_target_area);
}

template<unsigned DIM>
double DeathClockTargetAreaModifier<DIM>::GetAgeToStartGrowing()
{
    return mAgeToStartGrowing;
}

template<unsigned DIM>
void DeathClockTargetAreaModifier<DIM>::SetAgeToStartGrowing(double ageToStartGrowing)
{
    assert(ageToStartGrowing >= 0.0);
    mAgeToStartGrowing = ageToStartGrowing;
}

template<unsigned DIM>
double DeathClockTargetAreaModifier<DIM>::GetGrowthRate()
{
    return mGrowthRate;
}

template<unsigned DIM>
void DeathClockTargetAreaModifier<DIM>::SetGrowthRate(double growthRate)
{
    assert(growthRate >= 0.0);
    mGrowthRate = growthRate;
}

template<unsigned DIM>
void DeathClockTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AgeToStartGrowing>" << mAgeToStartGrowing << "</AgeToStartGrowing>\n";
    *rParamsFile << "\t\t\t<GrowthRate>" << mGrowthRate << "</GrowthRate>\n";

    // Next, call method on direct parent class
    AbstractTargetAreaModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DeathClockTargetAreaModifier<1>;
template class DeathClockTargetAreaModifier<2>;
template class DeathClockTargetAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeathClockTargetAreaModifier)
