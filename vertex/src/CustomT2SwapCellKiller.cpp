#include "CustomT2SwapCellKiller.hpp"

template<unsigned DIM>
CustomT2SwapCellKiller<DIM>::CustomT2SwapCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
    // Throw an exception if the population is not a CustomVertexBasedCellPopulation
    if (dynamic_cast<CustomVertexBasedCellPopulation<DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("A CustomT2SwapCellKiller should only be used together with a CustomVertexBasedCellPopulation.");
    }
}

template<unsigned DIM>
void CustomT2SwapCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    /*
     * This killer is different to other killers: it does not only check and label
     * cells for apoptosis or death, it actually carries out vertex rearrangements
     * and removes elements from the vertex mesh.
     *
     * We start with carrying out T2 swaps. Get the mesh and an element map.
     * The static_cast will work since we already know it's a CustomVertexBasedCellPopulation.
     */
    MutableVertexMesh<DIM,DIM>& mesh = static_cast<MutableVertexMesh<DIM,DIM>&>(this->mpCellPopulation->rGetMesh());
    CustomVertexBasedCellPopulation<DIM>* p_vertex_population = static_cast<CustomVertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);
    VertexElementMap element_map(mesh.GetNumAllElements());

    bool recheck_mesh = true;
    while (recheck_mesh == true)
    {
        // Note that whenever we call CheckForT2Swaps(), the element indices must run from zero up to mElements.size()-1
        recheck_mesh = mesh.CheckForT2Swaps(element_map);
        /*
         * There might have maximally one T2 swap happened above, where a vertex element was removed from the
         * mesh but the associated cell is still there. Here we check whether a new cell
         * underwent a T2 swap and label it as dead as well as record its location and ID.
         */
        for (unsigned elem_index = 0; elem_index < element_map.Size(); elem_index++)
        {
            CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(elem_index);
            if (element_map.IsDeleted(elem_index) && !(p_cell->IsDead()))
            {
                p_vertex_population->AddLocationOfT2Swap(mesh.GetLastT2SwapLocation());
                p_vertex_population->AddCellIdOfT2Swap(p_cell->GetCellId());
                p_vertex_population->AddAncestorOfT2Swap(p_cell->GetAncestor());
                p_vertex_population->AddIsApoptoticOfT2Swap(p_cell->HasApoptosisBegun());
                p_cell->Kill();

                // There can't have been more than one new cell death, so leave the for loop here.
                break;
            }
        }
    }

}

template<unsigned DIM>
void CustomT2SwapCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class CustomT2SwapCellKiller<1>;
template class CustomT2SwapCellKiller<2>;
template class CustomT2SwapCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CustomT2SwapCellKiller)
