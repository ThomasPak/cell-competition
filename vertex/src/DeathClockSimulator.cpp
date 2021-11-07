#include "SmartPointers.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"

#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "CustomUniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellAncestor.hpp"

#include "NagaiHondaForce.hpp"
#include "FarhadifarForce.hpp"

#include "CellProliferativePhasesWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"

#include "LogFile.hpp"

#include "RandomNumberGenerator.hpp"
#include "ExecutableSupport.hpp"
#include "SingletonSupport.hpp"
#include "OutputFileHandler.hpp"

#include "OffLatticeSimulationWithMinimumMaximumCellCount.hpp"
#include "DeathClockTargetAreaModifier.hpp"
#include "DeathClockCellKiller.hpp"
#include "DeathClockModifier.hpp"
#include "DeathClockSrnModel.hpp"

#include "PerimeterRatioWriter.hpp"
#include "T2SwapAncestorIsApoptoticWriter.hpp"
#include "CustomVertexBasedCellPopulation.hpp"

#include "DeathClockParameters.hpp"
#include "DeathClockSimulator.hpp"

int DeathClockSimulator(int argc, char *argv[],
        const std::string& input_string, std::string& output_string)
{
    // You should put all the main code within a try-catch, to ensure that
    // you clean up PETSc before quitting.
    try
    {
        // Process arguments
        if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help" ))
        {
            std::cerr << "Usage: " << argv[0] << " [KEY=VALUE]...\n";
            std::cerr << std::endl;
            std::cerr << DeathClockParameters::Help() << std::endl;
            return ExecutableSupport::EXIT_OK;
        }

        std::string parameter_string;
        for (int i = 1; i < argc; i++)
        {
            parameter_string += argv[i];
            parameter_string += '\n';
        }

        // Process input string
        parameter_string += input_string;
        Dictionary dict = create_dictionary(parameter_string);

        // Important: constructor of DeathClockParameters uses
        // RandomNumberGenerator to generate random birth times, so seed must
        // be read from the input before DeathClockParameters is constructed.
        unsigned seed = 0;
        if (dict.count("seed") == 1)
            seed = std::stoul(dict.at("seed"));

        SingletonSupport::SetupSingletons(seed);

        // Create DeathClockParameters object
        DeathClockParameters param;
        param.Update(dict);

        // Setup LogFile
        std::string output_path =
            param.m_output_directory + "/" + param.m_simulation_id + "/";

        LogFile* p_log = LogFile::Instance();
        p_log->Set(param.m_log_level, output_path);
        p_log->WriteHeader(output_path);

        // Print output
        std::cerr << param.Print() << std::endl;

        // Save output to file
        OutputFileHandler output_file_handler(output_path, false);

        out_stream death_clock_parameters_file =
            output_file_handler.OpenOutputFile("death_clock.parameters");
        *death_clock_parameters_file << param.Print();
        death_clock_parameters_file->close();

        // Run simulation
        RunDeathClockSimulation(param);

        // Destroy singletons
        SingletonSupport::DestroySingletons();

        // Empty output string
        output_string.assign("");

        return ExecutableSupport::EXIT_OK;
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());

        // Ensure singletons are destroyed
        SingletonSupport::DestroySingletons();

        // Write result
        output_string.assign("");

        return ExecutableSupport::EXIT_ERROR;
    }
}

void RunDeathClockSimulation(const DeathClockParameters& param)
{
    // Check whether param is valid
    if (!param.IsValid())
    {
        EXCEPTION("DeathClockParameters object is not valid");
    }

    // Use the honeycomb vertex mesh generator to create a vertex mesh.
    HoneycombVertexMeshGenerator generator(param.m_num_cells_across, param.m_num_cells_up);
    MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
    p_mesh->SetCellRearrangementThreshold(param.m_cell_rearrangement_threshold);
    p_mesh->SetCellRearrangementRatio(param.m_cell_rearrangement_ratio);
    p_mesh->SetT2Threshold(param.m_t2_threshold);
    p_mesh->SetCheckForInternalIntersections(true);

    // Define mutation state and proliferative type
    MAKE_PTR(WildTypeCellMutationState, p_state);
    MAKE_PTR(TransitCellProliferativeType, p_transit_type);

    // Create cells
    std::vector<CellPtr> cells;
    for (unsigned i = 0; i < param.m_initial_num_cells; i++)
    {
        // Get individual parameters
        CellCycleModel ccm = static_cast<CellCycleModel>(param.m_cell_cycle_models[i]);
        double g1_duration = param.m_g1_durations[i];
        double g2_duration = param.m_g2_durations[i];
        double uniform_g1_range = param.m_uniform_g1_ranges[i];
        unsigned max_transit_generation = param.m_max_transit_generations[i];
        double initial_tau = param.m_initial_taus[i];
        unsigned ancestor = param.m_ancestors[i];
        double reference_target_area = param.m_reference_target_areas[i];
        double apoptosis_time = param.m_apoptosis_times[i];
        double birth_time = param.m_birth_times[i];

        // Create cell cycle model
        AbstractSimpleGenerationalCellCycleModel* p_cell_cycle_model = nullptr;
        switch (ccm)
        {
            case exponential:
                p_cell_cycle_model = new ExponentialG1GenerationalCellCycleModel();
                break;

            case uniform:
                p_cell_cycle_model = new CustomUniformG1GenerationalCellCycleModel();
                break;

            default:
                NEVER_REACHED;
        }

        // Set cell cycle phase durations
        p_cell_cycle_model->SetMDuration(1e-12);
        p_cell_cycle_model->SetSDuration(1e-12);

        p_cell_cycle_model->SetTransitCellG1Duration(g1_duration);
        p_cell_cycle_model->SetG2Duration(g2_duration);

        if (ccm == uniform)
        {
            CustomUniformG1GenerationalCellCycleModel* p_uniform_ccm =
                dynamic_cast<CustomUniformG1GenerationalCellCycleModel*>(p_cell_cycle_model);
            assert(p_uniform_ccm);
            p_uniform_ccm->SetRange(uniform_g1_range);
        }

        // Set max transit generations
        p_cell_cycle_model->SetMaxTransitGenerations(max_transit_generation);

        // Create SRN model
        DeathClockSrnModel* p_srn_model = new DeathClockSrnModel;

        // Set initial conditions of SRN model
        std::vector<double> initial_conditions;
        initial_conditions.push_back(initial_tau);
        p_srn_model->SetInitialConditions(initial_conditions);

        // Create cell
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));

        // Set cell proliferative type
        p_cell->SetCellProliferativeType(p_transit_type);

        // Set cell ancestor
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (ancestor));
        p_cell->SetAncestor(p_cell_ancestor);

        // Set reference target area
        p_cell->GetCellData()->SetItem("reference target area", reference_target_area);

        // Set cell apoptosis time
        p_cell->SetApoptosisTime(apoptosis_time);

        // Set cell birth time
        p_cell->SetBirthTime(birth_time);

        // Push cell
        cells.push_back(p_cell);
    }

    // Define cell population
    CustomVertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

    // Convert map to unordered_map
    std::unordered_map<unsigned, unsigned> min_cell_count_for_ancestor;
    std::unordered_map<unsigned, unsigned> max_cell_count_for_ancestor;
    for (const auto& item : param.m_min_cell_count_for_ancestor)
    {
        min_cell_count_for_ancestor[item.first] = item.second;
    }
    for (const auto& item : param.m_max_cell_count_for_ancestor)
    {
        max_cell_count_for_ancestor[item.first] = item.second;
    }

    // Define simulation
    OffLatticeSimulationWithMinimumMaximumCellCount simulator(cell_population,
            param.m_min_cell_count,
            param.m_max_cell_count,
            min_cell_count_for_ancestor,
            max_cell_count_for_ancestor);

    // mpOdeSystem is initialised in constructor of OffLatticeSimulation,
    // so assign death clock parameters after construction.
    unsigned i = 0;
    for (auto p_cell = cell_population.Begin();
            p_cell != cell_population.End();
            ++p_cell)
    {
        // Get individual parameters
        double death_threshold = param.m_death_thresholds[i];

        double base_rate = param.m_base_rates[i];

        double neighbours_constant = param.m_neighbours_constants[i];
        double neighbours_g2_constant = param.m_neighbours_g2_constants[i];
        double neighbours_apoptosis_constant =
            param.m_neighbours_apoptosis_constants[i];
        double normalised_neighbours_g2_constant =
            param.m_normalised_neighbours_g2_constants[i];
        double normalised_neighbours_apoptosis_constant =
            param.m_normalised_neighbours_apoptosis_constants[i];

        double shared_edge_constant = param.m_shared_edge_constants[i];
        double shared_edge_g2_constant = param.m_shared_edge_g2_constants[i];
        double shared_edge_apoptosis_constant =
            param.m_shared_edge_apoptosis_constants[i];
        double normalised_shared_edge_g2_constant =
            param.m_normalised_shared_edge_g2_constants[i];
        double normalised_shared_edge_apoptosis_constant =
            param.m_normalised_shared_edge_apoptosis_constants[i];

        // Set cell death threshold
        p_cell->GetCellData()->SetItem("death threshold", death_threshold);

        // Get SRN model
        DeathClockSrnModel* p_srn_model =
            dynamic_cast<DeathClockSrnModel*>(p_cell->GetSrnModel());
        assert(p_srn_model);

        // Set parameters of SRN model
        p_srn_model->SetBaseRate(base_rate);

        p_srn_model->SetNeighboursConstant(neighbours_constant);
        p_srn_model->SetNeighboursG2Constant(neighbours_g2_constant);
        p_srn_model->SetNeighboursApoptosisConstant(neighbours_apoptosis_constant);
        p_srn_model->SetNormalisedNeighboursG2Constant(normalised_neighbours_g2_constant);
        p_srn_model->SetNormalisedNeighboursApoptosisConstant(normalised_neighbours_apoptosis_constant);

        p_srn_model->SetSharedEdgeConstant(shared_edge_constant);
        p_srn_model->SetSharedEdgeG2Constant(shared_edge_g2_constant);
        p_srn_model->SetSharedEdgeApoptosisConstant(shared_edge_apoptosis_constant);
        p_srn_model->SetNormalisedSharedEdgeG2Constant(normalised_shared_edge_g2_constant);
        p_srn_model->SetNormalisedSharedEdgeApoptosisConstant(normalised_shared_edge_apoptosis_constant);

        // Increment index
        i++;
    }

    // Add cell writers
    cell_population.AddCellWriter<CellProliferativePhasesWriter>();
    cell_population.AddCellWriter<CellAgesWriter>();
    cell_population.AddCellWriter<CellAncestorWriter>();

    // Add population writer
    cell_population.AddPopulationWriter<PerimeterRatioWriter>();
    cell_population.AddPopulationWriter<T2SwapAncestorIsApoptoticWriter>();

    // Configure simulator
    std::string output_directory = param.m_output_directory + '/' +
        param.m_simulation_id;
    simulator.SetOutputDirectory(output_directory);

    // Add force
    boost::shared_ptr<AbstractForce<2>> p_force;
    if (param.m_force_model == "farhadifar")
    {
        p_force.reset(new FarhadifarForce<2>);
    }
    else if (param.m_force_model == "nagai-honda")
    {
        p_force.reset(new NagaiHondaForce<2>);
    }
    else
    {
        std::string error_msg;
        error_msg += "Invalid force model (\"";
        error_msg += param.m_force_model;
        error_msg += "\")";
        EXCEPTION(error_msg);
    }

    // Configure FarhadifarForce
    auto p_farhadifar_force =
        boost::dynamic_pointer_cast<FarhadifarForce<2>>(p_force);

    if (p_farhadifar_force)
    {
        p_farhadifar_force->SetBoundaryLineTensionParameter(param.m_boundary_tension_parameter);
    }

    // Add force to simulator
    simulator.AddForce(p_force);

    // Add modifiers
    MAKE_PTR(DeathClockModifier<2>, p_death_clock_modifier);
    MAKE_PTR(DeathClockTargetAreaModifier<2>, p_growth_modifier);

    simulator.AddSimulationModifier(p_death_clock_modifier);
    simulator.AddSimulationModifier(p_growth_modifier);

    // Add cell killers
    MAKE_PTR_ARGS(DeathClockCellKiller<2>, p_death_clock_cell_killer,
            (&cell_population));
    simulator.AddCellKiller(p_death_clock_cell_killer);

    // Set simulation parameters
    simulator.SetEndTime(param.m_simulation_time);
    simulator.SetDt(param.m_dt);
    simulator.SetSamplingTimestepMultiple(param.m_sampling_timestep_multiple);

    // Solve
    std::cerr << "Solving...\n";
    simulator.Solve();
    std::cerr << "Done solving!\n";
}
