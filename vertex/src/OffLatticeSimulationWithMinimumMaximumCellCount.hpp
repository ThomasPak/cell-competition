#ifndef OFFLATTICESIMULATIONWITHMINIMUMMAXIMUMCELLCOUNT_HPP
#define OFFLATTICESIMULATIONWITHMINIMUMMAXIMUMCELLCOUNT_HPP

#include "OffLatticeSimulation.hpp"

#include <unordered_map>
#include <vector>

/**
 * Simple subclass of OffLatticeSimulation which just overloads
 * StoppingEventHasOccurred for stopping the simulation when the number of
 * cells hits a configurable minimum or maximum cell count.
 */
class OffLatticeSimulationWithMinimumMaximumCellCount : public OffLatticeSimulation<2>
{
private:

    /** Minimum cell count */
    unsigned mMinimumCellCount;

    /** Maximum cell count */
    unsigned mMaximumCellCount;

    /** Minimum cell count for ancestor */
    std::unordered_map<unsigned, unsigned> mMinimumCellCountForAncestor;

    /** Maximum cell count for ancestor */
    std::unordered_map<unsigned, unsigned> mMaximumCellCountForAncestor;

    /** Ancestors with minimum cell count */
    std::vector<unsigned> mAncestorsWithMinimumCellCount;

    /** Ancestors with maximum cell count */
    std::vector<unsigned> mAncestorsWithMaximumCellCount;

    /** Terminating */
    bool mTerminating = false;

    /** Define a stopping event which says stop if t>3.14 */
    bool StoppingEventHasOccurred() override;

public:
    OffLatticeSimulationWithMinimumMaximumCellCount (
            AbstractCellPopulation<2>& rCellPopulation,
            unsigned minimum_cell_count = 0u,
            unsigned maximum_cell_count = static_cast<unsigned>(-1),
            std::unordered_map<unsigned, unsigned> minimum_cell_count_for_ancestor = {},
            std::unordered_map<unsigned, unsigned> maximum_cell_count_for_ancestor = {}
            );

    /** Get minimum cell count */
    unsigned GetMinimumCellCount() const;

    /** Get maximum cell count */
    unsigned GetMaximumCellCount() const;

    /** Get ancestors with minimum cell count */
    std::vector<unsigned> GetAncestorsWithMinimumCellCount() const;

    /** Get minimum cell count for ancestor */
    unsigned GetMinimumCellCountForAncestor(unsigned ancestor) const;

    /** Get ancestors with maximum cell count */
    std::vector<unsigned> GetAncestorsWithMaximumCellCount() const;

    /** Get maximum cell count for ancestor */
    unsigned GetMaximumCellCountForAncestor(unsigned ancestor) const;
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(OffLatticeSimulationWithMinimumMaximumCellCount)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a OffLatticeSimulationWithMinimumMaximumCellCount.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulationWithMinimumMaximumCellCount * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* p_cell_population = &(t->rGetCellPopulation());
    const unsigned max_cell_count = t->GetMaximumCellCount();
    const unsigned min_cell_count = t->GetMinimumCellCount();

    ar & p_cell_population;
    ar & max_cell_count;
    ar & min_cell_count;
}

/**
 * De-serialize constructor parameters and initialise a OffLatticeSimulationWithMinimumMaximumCellCount.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, OffLatticeSimulationWithMinimumMaximumCellCount * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    unsigned max_cell_count;
    unsigned min_cell_count;

    ar >> p_cell_population;
    ar >> max_cell_count;
    ar >> min_cell_count;

    // Invoke inplace constructor to initialise instance
    ::new(t)OffLatticeSimulationWithMinimumMaximumCellCount(
            *p_cell_population,
            min_cell_count,
            max_cell_count);
}
}
} // namespace


#endif // OFFLATTICESIMULATIONWITHMINIMUMMAXIMUMCELLCOUNT_HPP
