#ifndef DEATHCLOCKCELLKILLER_HPP_
#define DEATHCLOCKCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A cell killer that kills cells if their death clock goes past the death threshold.
 */
template<unsigned DIM>
class DeathClockCellKiller : public AbstractCellKiller<DIM>
{
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

    public:

        /**
         * Default constructor.
         *
         * @param pCellPopulation pointer to the cell population
         */
        DeathClockCellKiller(AbstractCellPopulation<DIM>* pCellPopulation);

        /**
         * Loop over cells and start apoptosis if death clock has reached the
         * death threshold.
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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeathClockCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a DeathClockCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const DeathClockCellKiller<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a DeathClockCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, DeathClockCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)DeathClockCellKiller<DIM>(p_cell_population);
}
}
} // namespace ...

#endif // DEATHCLOCKCELLKILLER_HPP_
