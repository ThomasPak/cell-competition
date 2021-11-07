#ifndef CUSTOMUNIFORMG1GENERATIONALCELLCYCLEMODEL_HPP
#define CUSTOMUNIFORMG1GENERATIONALCELLCYCLEMODEL_HPP

#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"


/**
 * A cell cycle model where the G1 duration is drawn from an uniform random
 * distribution.  The range of this distribution is defined by the member
 * variable mRange.  The midpoint of this distribution is given by
 * GetStemCellG1Duration() and GetTransitCellG1Duration() for stem and transit
 * cells, respectively.
 *
 * The default value of mRange and is set to 2.
 *
 * It is required that the midpoint m satisfies m - 0.5 * mRange >= 0.
 */
class CustomUniformG1GenerationalCellCycleModel : public AbstractSimpleGenerationalCellCycleModel
{

private:

    /**
     * The range of the uniform distribution.
     * This parameter is initialised to 2.
     */
    double mRange;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationalCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
    }

protected:

    /**
     * Stochastically set the G1 duration following an uniform distribution.
     * Called on cell creation at the start of a simulation, and for both
     * parent and daughter cells at cell division.
     *
     * If the cell is differentiated, then the G1 phase duration is set to DBL_MAX,
     * so that the cell will never reach the end of G1 phase.
     */
    void SetG1Duration();

    /**
     * Protected copy-constructor for use by CreateCellCycleModel().
     *
     * The only way for external code to create a copy of a cell-cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    CustomUniformG1GenerationalCellCycleModel(const CustomUniformG1GenerationalCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    CustomUniformG1GenerationalCellCycleModel();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Return the range of the uniform distribution.
     *
     * @return mRange, the range of the distribution
     */
    double GetRange();

    /**
     * Set the range of the uniform distribution.
     *
     * @param range
     */
    void SetRange(double range);

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CustomUniformG1GenerationalCellCycleModel)


#endif // CUSTOMUNIFORMG1GENERATIONALCELLCYCLEMODEL_HPP
