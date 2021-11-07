#ifndef DEATHCLOCKSRNMODEL_HPP_
#define DEATHCLOCKSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "DeathClockOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

/**
 * A subclass of AbstractOdeSrnModel that includes a Death Clock ODE system in
 * the sub-cellular reaction network.
 */
class DeathClockSrnModel : public AbstractOdeSrnModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:
    /** Protected copy-constructor for use by CreateSrnModel().
     *
     * The only way for external code to create a copy of a SRN model is by
     * calling that method, to ensure that a model of the correct subclass is
     * created.  This copy-constructor helps subclasses to ensure that all
     * member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a
     * daughter cell upon cell division.  Note that the parent SRN model will
     * have had ResetForDivision() called just before CreateSrnModel() is
     * called, so performing an exact copy of the parent is suitable behaviour.
     * Any daughter-cell-specific initialisation can be done in
     * InitialiseDaughterCell().
     *
     * @param rModel  the SRN model to copy.
     */
    DeathClockSrnModel(const DeathClockSrnModel& rModel);

    /* The first public method is a constructor, which just calls the base
     * constructor.  Note you can include an optional argument to specify the ODE solver.*/
public:

    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver
     * object (allows the use of different ODE solvers)
     */
    DeathClockSrnModel(
            boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver
            = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of
     * this SRN model.
     *
     * @return a copy of the current SRN model.
     */
    AbstractSrnModel* CreateSrnModel();

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This overridden method sets up a new Death Clock ODE system.
     */
    void Initialise() override;

    /**
     * Overridden SimulateToCurrentTime() method for custom behaviour.
     *
     * This runs the ODE to the current time and then saves tau in the
     * CellData.
     */
    void SimulateToCurrentTime();

    /**
     * Output SRN model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSrnModelParameters(out_stream& rParamsFile);

    /**
     * Reset death clock before dividing.
     */
    void ResetForDivision() override;

    /**
     * Set tau.
     */
    void SetTau(double tau);

    /**
     * Set base rate.
     */
    void SetBaseRate(double base_rate);

    /**
     * Set neighbours constant.
     */
    void SetNeighboursConstant(double neighbours_constant);

    /**
     * Set neighbours G2 constant.
     */
    void SetNeighboursG2Constant(double neighbours_g2_constant);

    /**
     * Set neighbours apoptosis constant.
     */
    void SetNeighboursApoptosisConstant(double neighbours_apoptosis_constant);

    /**
     * Set normalised neighbours G2 constant.
     */
    void SetNormalisedNeighboursG2Constant(double normalised_neighbours_g2_constant);

    /**
     * Set normalised neighbours apoptosis constant.
     */
    void SetNormalisedNeighboursApoptosisConstant(double normalised_neighbours_apoptosis_constant);

    /**
     * Set shared edge constant.
     */
    void SetSharedEdgeConstant(double shared_edge_constant);

    /**
     * Set shared edge G2 constant.
     */
    void SetSharedEdgeG2Constant(double shared_edge_g2_constant);

    /**
     * Set shared edge apoptosis constant.
     */
    void SetSharedEdgeApoptosisConstant(double shared_edge_apoptosis_constant);

    /**
     * Set normalised shared edge G2 constant.
     */
    void SetNormalisedSharedEdgeG2Constant(double normalised_shared_edge_g2_constant);

    /**
     * Set normalised shared edge apoptosis constant.
     */
    void SetNormalisedSharedEdgeApoptosisConstant(double normalised_shared_edge_apoptosis_constant);

    /**
     * Update the inputs to the death clock ODE from CellData.
     */
    void UpdateDeathClockInputs();

    /**
     * @return the tau in this cell.
     */
    double GetTau() const;

    /**
     * @return the base rate of the cell.
     */
    double GetBaseRate() const;

    /**
     * @return the neighbours constant of the cell.
     */
    double GetNeighboursConstant() const;

    /**
     * @return the neighbours G2 constant of the cell.
     */
    double GetNeighboursG2Constant() const;

    /**
     * @return the neighbours apoptosis constant of the cell.
     */
    double GetNeighboursApoptosisConstant() const;

    /**
     * @return the normalised neighbours G2 constant of the cell.
     */
    double GetNormalisedNeighboursG2Constant() const;

    /**
     * @return the normalised neighbours apoptosis constant of the cell.
     */
    double GetNormalisedNeighboursApoptosisConstant() const;

    /**
     * @return the shared edge constant of the cell.
     */
    double GetSharedEdgeConstant() const;

    /**
     * @return the shared edge G2 constant of the cell.
     */
    double GetSharedEdgeG2Constant() const;

    /**
     * @return the shared edge apoptosis constant of the cell.
     */
    double GetSharedEdgeApoptosisConstant() const;

    /**
     * @return the normalised shared edge G2 constant of the cell.
     */
    double GetNormalisedSharedEdgeG2Constant() const;

    /**
     * @return the normalised shared edge apoptosis constant of the cell.
     */
    double GetNormalisedSharedEdgeApoptosisConstant() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DeathClockSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeathClockSrnModel)

#endif // DEATHCLOCKSRNMODEL_HPP_
