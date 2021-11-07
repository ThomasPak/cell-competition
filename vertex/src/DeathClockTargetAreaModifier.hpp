#ifndef DEATHCLOCKTARGETAREAMODIFIER_HPP_
#define DEATHCLOCKTARGETAREAMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractTargetAreaModifier.hpp"
#include "VertexBasedCellPopulation.hpp"

/**
 * This class modifies the cells' target areas in the same way as
 * TargetAreaLinearGrowthModifier, but assumes that the reference target area
 * is stored in each cell's CellData in the field "reference target area".
 *
 * /todo DeathClockTargetAreaModifier does not need to be derived from
 * AbstractTargetAreaModifier because it does not use the class member
 * "mReferenceTargetArea".
 */
template<unsigned DIM>
class DeathClockTargetAreaModifier : public AbstractTargetAreaModifier<DIM>
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
        archive & boost::serialization::base_object<AbstractTargetAreaModifier<DIM> >(*this);
        archive & mAgeToStartGrowing;
        archive & mGrowthRate;
    }

    /**
     * The age of a proliferating cell at which its target area should start growing.
     * Defaults to DOUBLE_UNSET. If this variable is set using SetAgeToStartGrowing(),
     * then it is used regardless of whether a phase-based cell-cycle model is used.
     */
    double mAgeToStartGrowing;

    /**
     * The growth rate of a proliferating cell's target area, when it is growing.
     * Defaults to DOUBLE_UNSET. This variable must be set if mAgeToStartGrowing is set,
     * otherwise an exception is thrown.
     */
    double mGrowthRate;

public:

    /**
     * Default constructor.
     */
    DeathClockTargetAreaModifier();

    /**
     * Destructor.
     */
    virtual ~DeathClockTargetAreaModifier() override;

    /**
     * Helper method to update the target area property of an individual cell.
     *
     * @param pCell pointer to a cell
     */
    void UpdateTargetAreaOfCell(const CellPtr pCell);

    /**
     * @return #mAgeToStartGrowing.
     */
    double GetAgeToStartGrowing();

    /**
     * Set #mAgeToStartGrowing.
     *
     * @param ageToStartGrowing the new value of #mAgeToStartGrowing
     */
    void SetAgeToStartGrowing(double ageToStartGrowing);

    /**
     * @return #mGrowthRate
     */
    double GetGrowthRate();

    /**
     * Set #mGrowthRate.
     *
     * @param growthRate the new value of #mGrowthRate
     */
    void SetGrowthRate(double growthRate);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeathClockTargetAreaModifier)

#endif // DEATHCLOCKTARGETAREAMODIFIER_HPP_
