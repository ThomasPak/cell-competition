#include <assert.h>
#include <unordered_map>
#include <vector>

#include "OffLatticeSimulationWithMinimumMaximumCellCount.hpp"

bool OffLatticeSimulationWithMinimumMaximumCellCount::StoppingEventHasOccurred()
{
    // Count cells by iterating over cell population
    unsigned cell_count = 0;

    std::unordered_map<unsigned, unsigned> cell_count_per_ancestor;
    for (auto ancestor : mAncestorsWithMinimumCellCount)
    {
        cell_count_per_ancestor[ancestor] = 0;
    }
    for (auto ancestor : mAncestorsWithMaximumCellCount)
    {
        cell_count_per_ancestor[ancestor] = 0;
    }

    for (auto p_cell = mrCellPopulation.Begin();
            p_cell != mrCellPopulation.End();
            ++p_cell)
    {
        // Increment cell count
        cell_count++;

        if (cell_count >= mMaximumCellCount)
        {
            mTerminating = true;
            break;
        }

        // Increment cell count per ancestor
        auto ancestor = p_cell->GetAncestor();
        if (cell_count_per_ancestor.count(ancestor) == 0)
        {
            cell_count_per_ancestor[ancestor] = 1;
        }
        else
        {
            cell_count_per_ancestor.at(ancestor)++;
        }

        for (auto ancestor : mAncestorsWithMaximumCellCount)
        {
            if (cell_count_per_ancestor[ancestor] >= mMaximumCellCountForAncestor[ancestor])
            {
                mTerminating = true;
                break;
            }
        }

        if (mTerminating == true)
        {
            break;
        }
    }

#ifndef NDEBUG
    unsigned cell_count_sum = 0;
    for (auto it : cell_count_per_ancestor)
    {
        cell_count_sum += it.second;
    }
    assert(cell_count == cell_count_sum);
#endif

    if (cell_count <= mMinimumCellCount)
    {
        mTerminating = true;
    }

    for (auto ancestor : mAncestorsWithMinimumCellCount)
    {
        if (cell_count_per_ancestor[ancestor] <= mMinimumCellCountForAncestor[ancestor])
        {
            mTerminating = true;
            break;
        }
    }

    // Check whether we are sampling results at the current timestep
    SimulationTime* p_time = SimulationTime::Instance();
    bool at_sampling_timestep = (p_time->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0);

    return mTerminating && at_sampling_timestep;
}

OffLatticeSimulationWithMinimumMaximumCellCount::OffLatticeSimulationWithMinimumMaximumCellCount(
        AbstractCellPopulation<2>& rCellPopulation,
        unsigned minimum_cell_count,
        unsigned maximum_cell_count,
        std::unordered_map<unsigned, unsigned> minimum_cell_count_for_ancestor,
        std::unordered_map<unsigned, unsigned> maximum_cell_count_for_ancestor
        )
    : OffLatticeSimulation<2>(rCellPopulation),
    mMinimumCellCount(minimum_cell_count),
    mMaximumCellCount(maximum_cell_count),
    mMinimumCellCountForAncestor(minimum_cell_count_for_ancestor),
    mMaximumCellCountForAncestor(maximum_cell_count_for_ancestor)
{
    for (auto it : mMinimumCellCountForAncestor)
    {
        mAncestorsWithMinimumCellCount.push_back(it.first);
    }

    for (auto it : mMaximumCellCountForAncestor)
    {
        mAncestorsWithMaximumCellCount.push_back(it.first);
    }
}

unsigned OffLatticeSimulationWithMinimumMaximumCellCount::GetMinimumCellCount() const
{
    return mMinimumCellCount;
}

unsigned OffLatticeSimulationWithMinimumMaximumCellCount::GetMaximumCellCount() const
{
    return mMaximumCellCount;
}

std::vector<unsigned>
OffLatticeSimulationWithMinimumMaximumCellCount::GetAncestorsWithMinimumCellCount()
    const
{
    return mAncestorsWithMinimumCellCount;
}

unsigned
OffLatticeSimulationWithMinimumMaximumCellCount::GetMinimumCellCountForAncestor(unsigned
        ancestor) const
{
    return mMinimumCellCountForAncestor.at(ancestor);
}

std::vector<unsigned>
OffLatticeSimulationWithMinimumMaximumCellCount::GetAncestorsWithMaximumCellCount()
    const
{
    return mAncestorsWithMaximumCellCount;
}

unsigned
OffLatticeSimulationWithMinimumMaximumCellCount::GetMaximumCellCountForAncestor(unsigned
        ancestor) const
{
    return mMaximumCellCountForAncestor.at(ancestor);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OffLatticeSimulationWithMinimumMaximumCellCount)
