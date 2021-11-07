#ifndef CUSTOMT2SWAPCELLKILLER_HPP
#define CUSTOMT2SWAPCELLKILLER_HPP

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "CustomVertexBasedCellPopulation.hpp"

/**
 * This killer functions identically to the T2SwapCellKiller, but in addition
 * the cell ancestor and death type for every T2 swap.  The death type is defined
 * as "extrusion" when the cell is not in apoptosis at the time of the T2 swap,
 * and as "apoptosis" when the cell is in apoptosis at the time of the T2 swap.
 */
template<unsigned DIM>
class CustomT2SwapCellKiller : public AbstractCellKiller<DIM>
{
    // The test will want to clean up the private member variable mpCellPopulation
    friend class TestCustomT2SwapCellKiller;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
    }

    /** Output file for death events. */
    out_stream mpDeathFile;

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
    CustomT2SwapCellKiller(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Overridden CheckAndLabelCellsForApoptosisOrDeath() method.
     *
     * Loop over cells and kill them if they are ready for a T2 swap.
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Overridden OutputCellKillerParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CustomT2SwapCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CustomT2SwapCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CustomT2SwapCellKiller<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a CustomT2SwapCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CustomT2SwapCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CustomT2SwapCellKiller<DIM>(p_cell_population);
}
} // namespace serialization
} // namespace boost

#endif // CUSTOMT2SWAPCELLKILLER_HPP
