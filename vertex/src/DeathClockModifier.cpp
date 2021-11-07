#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "VertexMesh.hpp"

#include "DeathClockModifier.hpp"

template<unsigned DIM>
DeathClockModifier<DIM>::DeathClockModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
DeathClockModifier<DIM>::~DeathClockModifier()
{
}

template<unsigned DIM>
void DeathClockModifier<DIM>::UpdateAtEndOfTimeStep(
        AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeathClockModifier<DIM>::SetupSolve(
        AbstractCellPopulation<DIM,DIM>& rCellPopulation,
        std::string outputDirectory)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DeathClockModifier<DIM>::OutputSimulationModifierParameters(
        out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned DIM>
void DeathClockModifier<DIM>::UpdateCellData(
        AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Ensure that mesh is vertex-based
    VertexMesh<DIM,DIM> *p_vm =
        dynamic_cast<VertexMesh<DIM,DIM>*>(&(rCellPopulation.rGetMesh()));

    if (p_vm == nullptr)
    {
        EXCEPTION("Mesh should be vertex-based");
    }

    // Next iterate over the population to compute and store each
    // cell's number of neighbours in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter =
            rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get my location index
        unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices =
            rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

        // Save number of neighbours and total edge length
        cell_iter->GetCellData()->SetItem("number of neighbours",
                static_cast<double>(neighbour_indices.size()));

        cell_iter->GetCellData()->SetItem("total edge length",
                p_vm->GetSurfaceAreaOfElement(elem_index));

        // Iterate over neighbours
        unsigned number_of_neighbours_apoptosis = 0;
        unsigned number_of_neighbours_g2 = 0;
        double shared_edge_length_apoptosis = 0.0;
        double shared_edge_length_g2 = 0.0;

        // Initialise shared_edge_length (this is not the same as total edge
        // length for boundary cells)
        double shared_edge_length = 0.0;

        for (auto iter = neighbour_indices.begin();
                iter != neighbour_indices.end();
                ++iter)
        {
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

            shared_edge_length += p_vm->GetEdgeLength(elem_index, *iter);

            if (p_cell->HasApoptosisBegun())
            {
                ++number_of_neighbours_apoptosis;
                shared_edge_length_apoptosis +=
                    p_vm->GetEdgeLength(elem_index, *iter);
            }

            // If cell cycle model is phase-based, record cells in G2
            if (AbstractPhaseBasedCellCycleModel *p_ccm =
                    dynamic_cast<AbstractPhaseBasedCellCycleModel*>(p_cell->GetCellCycleModel()))
            {
                if (p_ccm->GetCurrentCellCyclePhase() == G_TWO_PHASE)
                {
                    ++number_of_neighbours_g2;
                    shared_edge_length_g2 +=
                        p_vm->GetEdgeLength(elem_index, *iter);
                }
            }
        }

        cell_iter->GetCellData()->SetItem("number of neighbours apoptosis",
                static_cast<double>(number_of_neighbours_apoptosis));
        cell_iter->GetCellData()->SetItem("number of neighbours g2",
                static_cast<double>(number_of_neighbours_g2));
        cell_iter->GetCellData()->SetItem("shared edge length apoptosis",
                shared_edge_length_apoptosis);
        cell_iter->GetCellData()->SetItem("shared edge length g2",
                shared_edge_length_g2);
        cell_iter->GetCellData()->SetItem("shared edge length",
                shared_edge_length);
    }
}

// Explicit instantiation
template class DeathClockModifier<1>;
template class DeathClockModifier<2>;
template class DeathClockModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeathClockModifier)
