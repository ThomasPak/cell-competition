// This cell competition is adapted from Osborne 2017
// The paper tutorial for this paper can be found here:
// https://chaste.cs.ox.ac.uk/trac/wiki/PaperTutorials/CellBasedComparison2017

// I/O
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "FarhadifarDifferentialForce.hpp"
#include "CellAncestor.hpp"

// Project-specific force
#include "RandomMotionForce.hpp"

// Writers
#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "T2SwapAncestorIsApoptoticWriter.hpp"
#include "CustomVertexBasedCellPopulation.hpp"

#include "ExecutableSupport.hpp"
#include "SingletonSupport.hpp"
#include "CellLabelling.hpp"
#include "CellCompetitionParameters.hpp"
#include "CellCompetition.hpp"

int CellCompetitionSimulator(int argc, char *argv[],
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
            std::cerr << CellCompetitionParameters::Help() << std::endl;
            return ExecutableSupport::EXIT_OK;
        }

        std::string parameter_string;
        for (int i = 1; i < argc; i++)
        {
            parameter_string += argv[i];
            parameter_string += '\n';
        }

        // Process input string and create CellCompetitionParameters object
        parameter_string += input_string;
        Dictionary dict = create_dictionary(parameter_string);
        CellCompetitionParameters param;
        param.Update(dict);

        // Print output
        std::cerr << param.Print() << std::endl;

        // Setup singletons
        SingletonSupport::SetupSingletons(param.m_seed);

        // Run simulation
        RunCellCompetitionSimulation(param);

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

void RunCellCompetitionSimulation(const CellCompetitionParameters& param)
{
    // Create a simple 2D MutableVertexMesh
    HoneycombVertexMeshGenerator generator(param.m_num_cells_across,
            param.m_num_cells_up);
    MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
    p_mesh->SetCellRearrangementThreshold(param.m_cell_rearrangement_threshold);
    p_mesh->SetCellRearrangementRatio(param.m_cell_rearrangement_ratio);
    p_mesh->SetT2Threshold(param.m_t2_threshold);

    // Slows things down but can use a larger timestep and diffusion forces
    //p_mesh->SetCheckForInternalIntersections(true);

    // Set up cells, one for each VertexElement and with randomly generated
    // birth times sampled uniformly from [0, average cell cycle time]
    std::vector<CellPtr> cells;
    boost::shared_ptr<AbstractCellProperty>
        p_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
    CellsGenerator<ExponentialG1GenerationalCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

    for (unsigned i=0; i<cells.size(); i++)
    {
        // Set a target area rather than setting a growth modifier. (the
        // modifiers don't work correctly as making very long G1 phases)
        cells[i]->GetCellData()->SetItem("target area",
                param.m_cell_a_target_area);
    }

    // Create cell population
    CustomVertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

    // adjust cell cycle times
    for (auto cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
    {
        // Set ancestor
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (0)); // cell type A is ancestor 0
        cell_iter->SetAncestor(p_cell_ancestor);

        ExponentialG1GenerationalCellCycleModel* p_cell_cycle_model =
            static_cast<ExponentialG1GenerationalCellCycleModel*>(
                    cell_iter->GetCellCycleModel());

        // Set G1 duration (exponentially distributed)
        p_cell_cycle_model->SetTransitCellG1Duration(param.m_cell_a_g1_duration);

        // set M and S duration effectively to 0
        p_cell_cycle_model->SetMDuration(1e-12);
        p_cell_cycle_model->SetSDuration(1e-12);

        // Set G2 duration
        p_cell_cycle_model->SetG2Duration(param.m_cell_a_g2_duration);

        // Set maximum number of generations
        p_cell_cycle_model->SetMaxTransitGenerations(
                param.m_cell_a_max_generations);

        // Set birth time
        double birth_time = -
            p_cell_cycle_model->GetAverageTransitCellCycleTime() *
            RandomNumberGenerator::Instance()->ranf();
        (*cell_iter)->SetBirthTime(birth_time);

        p_cell_cycle_model->Initialise();
    }

    // Set population to output all data to results files
    cell_population.AddCellWriter<CellIdWriter>();
    cell_population.AddCellWriter<CellMutationStatesWriter>();
    cell_population.AddCellWriter<CellAncestorWriter>();

    //cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
    //cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();
    cell_population.AddPopulationWriter<T2SwapAncestorIsApoptoticWriter>();

    // We are using ParaView for visualization so we do not need the Chaste
    // visualizer data
    cell_population.SetOutputResultsForChasteVisualizer(false);

    // Set up cell-based simulation and output directory
    OffLatticeSimulation<2> simulator(cell_population);

    // Set up output directory
    std::string output_directory = param.m_output_directory + "/sim_" +
        param.m_simulation_id;
    simulator.SetOutputDirectory(output_directory);

    // Set time step and end time for simulation
    simulator.SetDt(param.m_dt);
    simulator.SetSamplingTimestepMultiple(param.m_sampling_timestep_multiple);
    simulator.SetEndTime(param.m_simulation_time);

    // Set up force law and pass it to the simulation
    MAKE_PTR(FarhadifarDifferentialForce<2>, p_force);

    p_force->SetCellAreaElasticityParameter(
            param.m_cell_a_elasticity_parameter);
    p_force->SetLabelledCellAreaElasticityParameter(
            param.m_cell_b_elasticity_parameter);

    p_force->SetCellPerimeterContractilityParameter(
            param.m_cell_a_contractility_parameter);
    p_force->SetLabelledCellPerimeterContractilityParameter(
            param.m_cell_b_contractility_parameter);


    p_force->SetCellCellLineTensionParameter(
            param.m_cell_a_line_tension_parameter);
    p_force->SetLabelledCellLabelledCellLineTensionParameter(
            param.m_cell_b_line_tension_parameter);

    p_force->SetLabelledCellCellLineTensionParameter(
            param.m_cell_ab_line_tension_parameter);

    p_force->SetCellBoundaryLineTensionParameter(
            param.m_cell_a_boundary_tension_parameter);
    p_force->SetLabelledCellBoundaryLineTensionParameter(
            param.m_cell_b_boundary_tension_parameter);

    simulator.AddForce(p_force);

    // Add some noise to avoid local minimum
    // Do not add force if parameter is zero
    if (param.m_random_movement_parameter > 0)
    {
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(param.m_random_movement_parameter);
        simulator.AddForce(p_random_force);
    }

    // Now label some cells
    boost::shared_ptr<AbstractCellProperty>
        p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
    if (param.m_random_labelling)
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state,
                param.m_initial_b_ratio, param.m_cell_b_target_area,
                param.m_cell_b_g1_duration,
                param.m_cell_b_g2_duration,
                param.m_cell_b_max_generations);
    else
        DeterministicallyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state,
                param.m_initial_b_ratio, param.m_cell_b_target_area,
                param.m_cell_b_g1_duration,
                param.m_cell_b_g2_duration,
                param.m_cell_b_max_generations);

    // Run simulation
    std::cerr << "Solving...\n";
    simulator.Solve();
    std::cerr << "Done solving!\n";

    // Count number of cells
    unsigned total_num_cells = cell_population.GetNumAllCells();
    unsigned num_labelled_cells = 0;
    for (auto elem_iter = cell_population.rGetMesh().GetElementIteratorBegin();
            elem_iter != cell_population.rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        CellPtr p_cell = cell_population.GetCellUsingLocationIndex(elem_index);

        if (p_cell->template HasCellProperty<CellLabel>())
            // This cell is labelled
            num_labelled_cells++;
    }

    double b_ratio =
        static_cast<double>(num_labelled_cells) /
        static_cast<double>(total_num_cells);

    // Open results.csv
    OutputFileHandler results_handler(output_directory, false);
    out_stream results_file = results_handler.OpenOutputFile("results.csv");

    // Output summary statistics to results file
    (*results_file) << "total-num-cells" << ","
                    << "b-ratio" << std::endl;

    (*results_file) << total_num_cells << ","
                    << b_ratio
                    << std::endl;

    // Tidy up
    results_file->close();
}
