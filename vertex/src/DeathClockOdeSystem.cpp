#include "DeathClockOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

DeathClockOdeSystem::DeathClockOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(1)
{
    mpSystemInfo = OdeSystemInformation<DeathClockOdeSystem>::Instance();

    SetDefaultInitialCondition(0, 0.0); // soon overwritten

    // Default parameter values are 0
    this->mParameters.push_back(0.0); // "number of neighbours"
    this->mParameters.push_back(0.0); // "number of neighbours g2"
    this->mParameters.push_back(0.0); // "number of neighbours apoptosis"

    this->mParameters.push_back(0.0); // "shared edge length"
    this->mParameters.push_back(0.0); // "shared edge length g2"
    this->mParameters.push_back(0.0); // "shared edge length apoptosis"

    this->mParameters.push_back(0.0); // "base rate"

    this->mParameters.push_back(0.0); // "neighbours constant"
    this->mParameters.push_back(0.0); // "neighbours g2 constant"
    this->mParameters.push_back(0.0); // "neighbours apoptosis constant"
    this->mParameters.push_back(0.0); // "normalised neighbours g2 constant"
    this->mParameters.push_back(0.0); // "normalised neighbours apoptosis constant"

    this->mParameters.push_back(0.0); // "shared edge constant"
    this->mParameters.push_back(0.0); // "shared edge g2 constant"
    this->mParameters.push_back(0.0); // "shared edge apoptosis constant"
    this->mParameters.push_back(0.0); // "normalised shared edge g2 constant"
    this->mParameters.push_back(0.0); // "normalised shared edge apoptosis constant"

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

DeathClockOdeSystem::~DeathClockOdeSystem()
{
}

void DeathClockOdeSystem::EvaluateYDerivatives(double time,
        const std::vector<double>& rY, std::vector<double>& rDY)
{
    // Get inputs
    double number_of_neighbours = this->GetParameter("number of neighbours");
    double number_of_neighbours_g2 = this->GetParameter("number of neighbours g2");
    double number_of_neighbours_apoptosis =
        this->GetParameter("number of neighbours apoptosis");

    double shared_edge_length = this->GetParameter("shared edge length");
    double shared_edge_length_g2 = this->GetParameter("shared edge length g2");
    double shared_edge_length_apoptosis =
        this->GetParameter("shared edge length apoptosis");

    // Get parameters
    double base_rate = this->GetParameter("base rate");

    double neighbours_constant = this->GetParameter("neighbours constant");
    double neighbours_g2_constant = this->GetParameter("neighbours g2 constant");
    double neighbours_apoptosis_constant =
        this->GetParameter("neighbours apoptosis constant");
    double normalised_neighbours_g2_constant = this->GetParameter("normalised neighbours g2 constant");
    double normalised_neighbours_apoptosis_constant =
        this->GetParameter("normalised neighbours apoptosis constant");

    double shared_edge_constant = this->GetParameter("shared edge constant");
    double shared_edge_g2_constant = this->GetParameter("shared edge g2 constant");
    double shared_edge_apoptosis_constant =
        this->GetParameter("shared edge apoptosis constant");
    double normalised_shared_edge_g2_constant = this->GetParameter("normalised shared edge g2 constant");
    double normalised_shared_edge_apoptosis_constant =
        this->GetParameter("normalised shared edge apoptosis constant");

    // Compute derivative
    rDY[0] = base_rate
        + neighbours_constant * number_of_neighbours
        + neighbours_g2_constant * number_of_neighbours_g2
        + neighbours_apoptosis_constant * number_of_neighbours_apoptosis
        + shared_edge_constant * shared_edge_length
        + shared_edge_g2_constant * shared_edge_length_g2
        + shared_edge_apoptosis_constant * shared_edge_length_apoptosis;


    if (number_of_neighbours > 0)
    {
        rDY[0] += normalised_neighbours_g2_constant * number_of_neighbours_g2 / number_of_neighbours
            + normalised_neighbours_apoptosis_constant * number_of_neighbours_apoptosis / number_of_neighbours;
    }

    if (shared_edge_length > 0)
    {
        rDY[0] += normalised_shared_edge_g2_constant * shared_edge_length_g2 / shared_edge_length
            + normalised_shared_edge_apoptosis_constant * shared_edge_length_apoptosis / shared_edge_length;
    }
}

template<>
void OdeSystemInformation<DeathClockOdeSystem>::Initialise()
{
    // Tau
    this->mVariableNames.push_back("tau");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);

    // Death clock base rate
    this->mParameterNames.push_back("base rate");
    this->mParameterUnits.push_back("dimensionless");

    // Neighbours constant
    this->mParameterNames.push_back("neighbours constant");
    this->mParameterUnits.push_back("dimensionless");

    // Number of neighbours
    this->mParameterNames.push_back("number of neighbours");
    this->mParameterUnits.push_back("dimensionless");

    // Neighbours G2 constant
    this->mParameterNames.push_back("neighbours g2 constant");
    this->mParameterUnits.push_back("dimensionless");

    // Number of neighbourss G2
    this->mParameterNames.push_back("number of neighbours g2");
    this->mParameterUnits.push_back("dimensionless");

    // Neighbours apoptosis constant
    this->mParameterNames.push_back("neighbours apoptosis constant");
    this->mParameterUnits.push_back("dimensionless");

    // Number of neighbourss apoptosis
    this->mParameterNames.push_back("number of neighbours apoptosis");
    this->mParameterUnits.push_back("dimensionless");

    // Normalised neighbours G2 constant
    this->mParameterNames.push_back("normalised neighbours g2 constant");
    this->mParameterUnits.push_back("dimensionless");

    // Normalised neighbours apoptosis constant
    this->mParameterNames.push_back("normalised neighbours apoptosis constant");
    this->mParameterUnits.push_back("dimensionless");

    // Shared edge constant
    this->mParameterNames.push_back("shared edge constant");
    this->mParameterUnits.push_back("dimensionless");

    // Shared edge length
    this->mParameterNames.push_back("shared edge length");
    this->mParameterUnits.push_back("dimensionless");

    // Shared edge G2 constant
    this->mParameterNames.push_back("shared edge g2 constant");
    this->mParameterUnits.push_back("dimensionless");

    // Shared edge length G2
    this->mParameterNames.push_back("shared edge length g2");
    this->mParameterUnits.push_back("dimensionless");

    // Shared edge apoptosis constant
    this->mParameterNames.push_back("shared edge apoptosis constant");
    this->mParameterUnits.push_back("dimensionless");

    // Shared edge length apoptosis
    this->mParameterNames.push_back("shared edge length apoptosis");
    this->mParameterUnits.push_back("dimensionless");

    // Normalised shared edge G2 constant
    this->mParameterNames.push_back("normalised shared edge g2 constant");
    this->mParameterUnits.push_back("dimensionless");

    // Normalised shared edge apoptosis constant
    this->mParameterNames.push_back("normalised shared edge apoptosis constant");
    this->mParameterUnits.push_back("dimensionless");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeathClockOdeSystem)
