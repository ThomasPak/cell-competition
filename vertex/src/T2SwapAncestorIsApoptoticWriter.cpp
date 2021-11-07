#include "Warnings.hpp"

#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CustomVertexBasedCellPopulation.hpp"

#include "T2SwapAncestorIsApoptoticWriter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
T2SwapAncestorIsApoptoticWriter<ELEMENT_DIM, SPACE_DIM>::T2SwapAncestorIsApoptoticWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("T2AncestorsIsApoptotics.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void T2SwapAncestorIsApoptoticWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void T2SwapAncestorIsApoptoticWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void T2SwapAncestorIsApoptoticWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void T2SwapAncestorIsApoptoticWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void T2SwapAncestorIsApoptoticWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // If the dynamic cast to CustomVertexBasedCellPopulation is successful,
    // use the specialised function.
    if (auto p_cell_pop = dynamic_cast<CustomVertexBasedCellPopulation<SPACE_DIM>*>(pCellPopulation))
    {
        Visit(p_cell_pop);
    }
    else
    {
        WARN_ONCE_ONLY("Using T2SwapAncestorIsApoptoticWriter with VertexBasedCellPopulation: "
                "apoptosis and extrusion events are not being counted.");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void T2SwapAncestorIsApoptoticWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CustomVertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    std::vector<unsigned> t2_ancestors = pCellPopulation->GetAncestorsOfT2Swaps();
    std::vector<bool> t2_is_apoptotics = pCellPopulation->GetIsApoptoticsOfT2Swaps();

    *this->mpOutStream << t2_ancestors.size() << "\t";
    assert( t2_ancestors.size() == t2_is_apoptotics.size());

    for (unsigned index = 0;  index < t2_ancestors.size(); index++)
    {
        *this->mpOutStream << t2_ancestors[index] << "\t";
        if (t2_is_apoptotics[index])
        {
            *this->mpOutStream << "apoptosis" << "\t";
        }
        else
        {
            *this->mpOutStream << "extrusion" << "\t";
        }
    }

    pCellPopulation->ClearAncestorsAndIsApoptoticsOfT2Swaps();
}

// Explicit instantiation
template class T2SwapAncestorIsApoptoticWriter<1,1>;
template class T2SwapAncestorIsApoptoticWriter<1,2>;
template class T2SwapAncestorIsApoptoticWriter<2,2>;
template class T2SwapAncestorIsApoptoticWriter<1,3>;
template class T2SwapAncestorIsApoptoticWriter<2,3>;
template class T2SwapAncestorIsApoptoticWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(T2SwapAncestorIsApoptoticWriter)
