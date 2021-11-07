#include "PerimeterRatioWriter.hpp"

#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PerimeterRatioWriter<ELEMENT_DIM, SPACE_DIM>::PerimeterRatioWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("perimeterratio.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PerimeterRatioWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    // Make sure the cell population is updated
    pCellPopulation->Update();

    // Initialise helper containers
    std::vector<double> cell_areas;
    std::vector<double> perimeter_ratios;
    std::set<unsigned> seen_ancestors;

    // Iterate over cells
    for (auto cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        // Get ancestor
        unsigned ancestor = cell_iter->GetAncestor();

        // Check if ancestor was set
        if (ancestor == UNSIGNED_UNSET)
        {
            EXCEPTION("Cell ancestor was not set");
        }

        // If not yet seen, add to seen_ancestors and resize cell_areas,
        // perimeter_ratios if needed
        if (seen_ancestors.count(ancestor) == 0)
        {
            seen_ancestors.insert(ancestor);
            if (cell_areas.size() < (ancestor + 1))
            {
                cell_areas.resize(ancestor + 1);
                perimeter_ratios.resize(ancestor + 1);
            }
        }

        // Get cell area
        double cell_area = pCellPopulation->GetVolumeOfCell(*cell_iter);

        // Initialise perimeter ratio
        double perimeter_ratio = 0.0;

        // Get the set of neighbouring element indices
        unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        std::set<unsigned> neighbour_elem_indices = pCellPopulation->rGetMesh().GetNeighbouringElementIndices(elem_index);

        // Iterate over these neighbours
        for (auto neighbour_iter = neighbour_elem_indices.begin();
             neighbour_iter != neighbour_elem_indices.end();
             ++neighbour_iter)
        {
            unsigned neighbour_index = *neighbour_iter;

            // Get the ancestor of this neighbour
            unsigned neighbour_ancestor =
                pCellPopulation->GetCellUsingLocationIndex(neighbour_index)->GetAncestor();

            if (ancestor == neighbour_ancestor)
            {
                // Get the length of the edge shared with this neighbour
                double edge_length =
                    pCellPopulation->rGetMesh().GetEdgeLength(elem_index,
                            neighbour_index);
                perimeter_ratio += edge_length;
            }

        }

        // Calculate perimeter_ratio
        perimeter_ratio /= pCellPopulation->rGetMesh().GetSurfaceAreaOfElement(elem_index);

        // Record cell area and weighted perimeter ratio
        cell_areas[ancestor] += cell_area;
        perimeter_ratios[ancestor] += perimeter_ratio * cell_area;
    }

    // Calculate and output weighted perimeter ratios
    unsigned index = 0;
    for (auto ancestor_iter = seen_ancestors.begin();
            ancestor_iter != seen_ancestors.end();
            ++ancestor_iter)
    {
        double weighted_perimeter_ratio = perimeter_ratios[*ancestor_iter] /
            cell_areas[*ancestor_iter];

        *this->mpOutStream << *ancestor_iter << " " << weighted_perimeter_ratio << "\t";

        index++;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PerimeterRatioWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PerimeterRatioWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PerimeterRatioWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PerimeterRatioWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

// Explicit instantiation
template class PerimeterRatioWriter<1,1>;
template class PerimeterRatioWriter<1,2>;
template class PerimeterRatioWriter<2,2>;
template class PerimeterRatioWriter<1,3>;
template class PerimeterRatioWriter<2,3>;
template class PerimeterRatioWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PerimeterRatioWriter)
