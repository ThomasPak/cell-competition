#include "AbstractVertexBasedDivisionRule.hpp"
#include "AbstractCellBasedSimulation.hpp"
#include "CustomT2SwapCellKiller.hpp"
#include "SmartPointers.hpp"

#include "CustomVertexBasedCellPopulation.hpp"

template<unsigned DIM>
CustomVertexBasedCellPopulation<DIM>::CustomVertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh,
                          std::vector<CellPtr>& rCells,
                          bool deleteMesh,
                          bool validate,
                          const std::vector<unsigned> locationIndices)
    : VertexBasedCellPopulation<DIM>(rMesh, rCells, deleteMesh, validate, locationIndices)
{
}

template<unsigned DIM>
CustomVertexBasedCellPopulation<DIM>::CustomVertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh)
    : VertexBasedCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
CustomVertexBasedCellPopulation<DIM>::~CustomVertexBasedCellPopulation()
{
}

template<unsigned DIM>
std::vector<unsigned> CustomVertexBasedCellPopulation<DIM>::GetAncestorsOfT2Swaps()
{
    return mAncestorsOfT2Swaps;
}

template<unsigned DIM>
std::vector<bool> CustomVertexBasedCellPopulation<DIM>::GetIsApoptoticsOfT2Swaps()
{
    return mIsApoptoticsOfT2Swaps;
}

template<unsigned DIM>
void CustomVertexBasedCellPopulation<DIM>::AddAncestorOfT2Swap(unsigned ancestorOfT2Swap)
{
    mAncestorsOfT2Swaps.push_back(ancestorOfT2Swap);
}

template<unsigned DIM>
void CustomVertexBasedCellPopulation<DIM>::AddIsApoptoticOfT2Swap(bool isApoptoticOfT2Swap)
{
    mIsApoptoticsOfT2Swaps.push_back(isApoptoticOfT2Swap);
}

template<unsigned DIM>
void CustomVertexBasedCellPopulation<DIM>::ClearAncestorsAndIsApoptoticsOfT2Swaps()
{
    mAncestorsOfT2Swaps.clear();
    mIsApoptoticsOfT2Swaps.clear();
}

template<unsigned DIM>
void CustomVertexBasedCellPopulation<DIM>::SimulationSetupHook(AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    MAKE_PTR_ARGS(CustomT2SwapCellKiller<DIM>, p_t2_swap_cell_killer, (this));
    pSimulation->AddCellKiller(p_t2_swap_cell_killer);
}

// Explicit instantiation
template class CustomVertexBasedCellPopulation<1>;
template class CustomVertexBasedCellPopulation<2>;
template class CustomVertexBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CustomVertexBasedCellPopulation)
