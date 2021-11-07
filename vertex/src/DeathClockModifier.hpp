#ifndef DEATHCLOCKMODIFIER_HPP_
#define DEATHCLOCKMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class that stores several quantities related to neighbouring
 * cells in CellData.  This includes:
 * - number of neighbours
 * - number of neighbours in G2 phase
 * - number of neighbours in apoptosis
 * - shared edge length (with neighbours)
 * - shared edge length in G2 phase
 * - shared edge length in apoptosis
 *
 * To be used in conjunction with Death Clock cell cycle models.
 *
 * TODO: future additions include
 * - store concentration of morphogen
 */
template <unsigned DIM>
class DeathClockModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    DeathClockModifier();

    /**
     * Destructor.
     */
    virtual ~DeathClockModifier() override;

    /*
     * Next, we override the {{{UpdateAtEndOfTimeStep()}}} method, which
     * specifies what to do to the simulation at the end of each time step. In
     * this class, we simply do nothing.
     */

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.  In
     * this class, we call UpdateCellData().
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation) override;

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.  In this
     * class, we call UpdateCellData().
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
            std::string outputDirectory) override;

    /**
     * Helper method to store the neighbours info in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeathClockModifier)

#endif // DEATHCLOCKMODIFIER_HPP_
