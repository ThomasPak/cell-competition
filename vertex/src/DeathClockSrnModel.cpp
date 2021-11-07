#include "DeathClockSrnModel.hpp"

DeathClockSrnModel::DeathClockSrnModel(const DeathClockSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * These lines copy the ODE system.
     */
    assert(rModel.GetOdeSystem());
    SetOdeSystem(new DeathClockOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));

    // Copy parameters
    SetBaseRate(rModel.GetBaseRate());

    SetNeighboursConstant(rModel.GetNeighboursConstant());
    SetNeighboursG2Constant(rModel.GetNeighboursG2Constant());
    SetNeighboursApoptosisConstant(rModel.GetNeighboursApoptosisConstant());
    SetNormalisedNeighboursG2Constant(rModel.GetNormalisedNeighboursG2Constant());
    SetNormalisedNeighboursApoptosisConstant(rModel.GetNormalisedNeighboursApoptosisConstant());

    SetSharedEdgeConstant(rModel.GetSharedEdgeConstant());
    SetSharedEdgeG2Constant(rModel.GetSharedEdgeG2Constant());
    SetSharedEdgeApoptosisConstant(rModel.GetSharedEdgeApoptosisConstant());
    SetNormalisedSharedEdgeG2Constant(rModel.GetNormalisedSharedEdgeG2Constant());
    SetNormalisedSharedEdgeApoptosisConstant(rModel.GetNormalisedSharedEdgeApoptosisConstant());
}

DeathClockSrnModel::DeathClockSrnModel(
        boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(1,
            boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<DeathClockSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
#else
        mpOdeSolver = CellCycleModelOdeSolver<DeathClockSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.1);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

AbstractSrnModel* DeathClockSrnModel::CreateSrnModel()
{
    return new DeathClockSrnModel(*this);
}

void DeathClockSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new DeathClockOdeSystem);
}

void DeathClockSrnModel::SimulateToCurrentTime()
{
    // Update inputs from CellData to ODE system
    UpdateDeathClockInputs();

    // run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();

    /* this line outputs the ODE system variable to CellData. */
    mpCell->GetCellData()->SetItem("tau", mpOdeSystem->rGetStateVariables()[0]);
}

void DeathClockSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

void DeathClockSrnModel::ResetForDivision()
{
    AbstractOdeSrnModel::ResetForDivision();
    std::vector<double>& state_variables = GetOdeSystem()->rGetStateVariables();
    state_variables[0] = 0.0;
}

void DeathClockSrnModel::SetTau(double tau)
{
    assert(mpOdeSystem != nullptr);
    std::vector<double>& state_variables = GetOdeSystem()->rGetStateVariables();
    state_variables[0] = tau;
}

void DeathClockSrnModel::SetBaseRate(double base_rate)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("base rate", base_rate);
}

void DeathClockSrnModel::SetNeighboursConstant(double neighbours_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("neighbours constant", neighbours_constant);
}

void DeathClockSrnModel::SetNeighboursG2Constant(double neighbours_g2_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("neighbours g2 constant", neighbours_g2_constant);
}

void DeathClockSrnModel::SetNeighboursApoptosisConstant(double neighbours_apoptosis_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("neighbours apoptosis constant", neighbours_apoptosis_constant);
}

void DeathClockSrnModel::SetNormalisedNeighboursG2Constant(
        double normalised_neighbours_g2_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("normalised neighbours g2 constant",
            normalised_neighbours_g2_constant);
}

void DeathClockSrnModel::SetNormalisedNeighboursApoptosisConstant(
        double normalised_neighbours_apoptosis_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("normalised neighbours apoptosis constant",
            normalised_neighbours_apoptosis_constant);
}

void DeathClockSrnModel::SetSharedEdgeConstant(double shared_edge_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("shared edge constant", shared_edge_constant);
}

void DeathClockSrnModel::SetSharedEdgeG2Constant(double shared_edge_g2_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("shared edge g2 constant", shared_edge_g2_constant);
}

void DeathClockSrnModel::SetSharedEdgeApoptosisConstant(double shared_edge_apoptosis_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("shared edge apoptosis constant", shared_edge_apoptosis_constant);
}

void DeathClockSrnModel::SetNormalisedSharedEdgeG2Constant(
        double normalised_shared_edge_g2_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("normalised shared edge g2 constant",
            normalised_shared_edge_g2_constant);
}

void DeathClockSrnModel::SetNormalisedSharedEdgeApoptosisConstant(
        double normalised_shared_edge_apoptosis_constant)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->SetParameter("normalised shared edge apoptosis constant",
            normalised_shared_edge_apoptosis_constant);
}

void DeathClockSrnModel::UpdateDeathClockInputs()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double number_of_neighbours = mpCell->GetCellData()->GetItem("number of neighbours");
    double number_of_neighbours_g2 =
        mpCell->GetCellData()->GetItem("number of neighbours g2");
    double number_of_neighbours_apoptosis =
        mpCell->GetCellData()->GetItem("number of neighbours apoptosis");

    mpOdeSystem->SetParameter("number of neighbours", number_of_neighbours);
    mpOdeSystem->SetParameter("number of neighbours g2",
            number_of_neighbours_g2);
    mpOdeSystem->SetParameter("number of neighbours apoptosis",
            number_of_neighbours_apoptosis);

    double shared_edge_length = mpCell->GetCellData()->GetItem("shared edge length");
    double shared_edge_length_g2 =
        mpCell->GetCellData()->GetItem("shared edge length g2");
    double shared_edge_length_apoptosis =
        mpCell->GetCellData()->GetItem("shared edge length apoptosis");

    mpOdeSystem->SetParameter("shared edge length", shared_edge_length);
    mpOdeSystem->SetParameter("shared edge length g2",
            shared_edge_length_g2);
    mpOdeSystem->SetParameter("shared edge length apoptosis",
            shared_edge_length_apoptosis);
}

double DeathClockSrnModel::GetTau() const
{
    assert(mpOdeSystem != nullptr);
    double tau = mpOdeSystem->rGetStateVariables()[0];
    return tau;
}

double DeathClockSrnModel::GetBaseRate() const
{
    assert(mpOdeSystem != nullptr);
    double base_rate = mpOdeSystem->GetParameter("base rate");
    return base_rate;
}

double DeathClockSrnModel::GetNeighboursConstant() const
{
    assert(mpOdeSystem != nullptr);
    double neighbours_constant = mpOdeSystem->GetParameter("neighbours constant");
    return neighbours_constant;
}

double DeathClockSrnModel::GetNeighboursG2Constant() const
{
    assert(mpOdeSystem != nullptr);
    double neighbours_g2_constant =
        mpOdeSystem->GetParameter("neighbours g2 constant");
    return neighbours_g2_constant;
}

double DeathClockSrnModel::GetNeighboursApoptosisConstant() const
{
    assert(mpOdeSystem != nullptr);
    double neighbours_apoptosis_constant =
        mpOdeSystem->GetParameter("neighbours apoptosis constant");
    return neighbours_apoptosis_constant;
}

double DeathClockSrnModel::GetNormalisedNeighboursG2Constant() const
{
    assert(mpOdeSystem != nullptr);
    double normalised_neighbours_g2_constant =
        mpOdeSystem->GetParameter("normalised neighbours g2 constant");
    return normalised_neighbours_g2_constant;
}

double DeathClockSrnModel::GetNormalisedNeighboursApoptosisConstant() const
{
    assert(mpOdeSystem != nullptr);
    double normalised_neighbours_apoptosis_constant =
        mpOdeSystem->GetParameter("normalised neighbours apoptosis constant");
    return normalised_neighbours_apoptosis_constant;
}

double DeathClockSrnModel::GetSharedEdgeConstant() const
{
    assert(mpOdeSystem != nullptr);
    double shared_edge_constant = mpOdeSystem->GetParameter("shared edge constant");
    return shared_edge_constant;
}

double DeathClockSrnModel::GetSharedEdgeG2Constant() const
{
    assert(mpOdeSystem != nullptr);
    double shared_edge_g2_constant =
        mpOdeSystem->GetParameter("shared edge g2 constant");
    return shared_edge_g2_constant;
}

double DeathClockSrnModel::GetSharedEdgeApoptosisConstant() const
{
    assert(mpOdeSystem != nullptr);
    double shared_edge_apoptosis_constant =
        mpOdeSystem->GetParameter("shared edge apoptosis constant");
    return shared_edge_apoptosis_constant;
}

double DeathClockSrnModel::GetNormalisedSharedEdgeG2Constant() const
{
    assert(mpOdeSystem != nullptr);
    double normalised_shared_edge_g2_constant =
        mpOdeSystem->GetParameter("normalised shared edge g2 constant");
    return normalised_shared_edge_g2_constant;
}

double DeathClockSrnModel::GetNormalisedSharedEdgeApoptosisConstant() const
{
    assert(mpOdeSystem != nullptr);
    double normalised_shared_edge_apoptosis_constant =
        mpOdeSystem->GetParameter("normalised shared edge apoptosis constant");
    return normalised_shared_edge_apoptosis_constant;
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeathClockSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeathClockSrnModel)
