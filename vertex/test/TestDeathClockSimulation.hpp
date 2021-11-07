#ifndef TESTDEATHCLOCK_HPP_
#define TESTDEATHCLOCK_HPP_

#include <cxxtest/TestSuite.h>

#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellAgesWriter.hpp"

#include "DeathClockSrnModel.hpp"
#include "DeathClockCellKiller.hpp"
#include "DeathClockModifier.hpp"

#include "PetscSetupAndFinalize.hpp"

/**
 * \todo TODO implement test for normalised g2, normalised apoptosis and shared
 * edge
 */
class TestDeathClock : public AbstractCellBasedTestSuite
{

private:

    template<class CellCycleModel>
    std::vector<CellPtr> create_cells(
            double g1_duration,
            double g2_duration,
            std::vector<double> birth_times,
            std::vector<double> initial_taus = std::vector<double>())
    {
        // If initial_taus is not empty, it must have the same number of
        // elements as birth_times
        assert(initial_taus.empty() || (birth_times.size() == initial_taus.size()));

        // Define mutation state and proliferative type
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_stem_type);

        // Initialise cells vector and loop
        std::vector<CellPtr> cells;
        for (unsigned i = 0; i < birth_times.size(); i++)
        {
            // Create cell cycle model
            CellCycleModel* p_cell_cycle_model = new CellCycleModel();

            // Set cell cycle phase durations
            p_cell_cycle_model->SetMDuration(1e-12);
            p_cell_cycle_model->SetSDuration(1e-12);

            p_cell_cycle_model->SetTransitCellG1Duration(g1_duration);
            p_cell_cycle_model->SetG2Duration(g2_duration);

            // Create SRN model
            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel;

            // Set initial conditions of SRN model
            std::vector<double> initial_conditions;
            initial_conditions.push_back(
                    initial_taus.empty() ? 0.0 : initial_taus[i]);
            p_srn_model->SetInitialConditions(initial_conditions);

            // Create cell
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));

            // Set cell proliferative type
            p_cell->SetCellProliferativeType(p_stem_type);

            // Set cell birth time
            p_cell->SetBirthTime(birth_times[i]);

            // Push cell
            cells.push_back(p_cell);
        }

        return cells;
    }

    void set_death_clock_parameters(
            double death_threshold,
            double base_rate,
            double neighbours_constant,
            double neighbours_g2_constant,
            double neighbours_apoptosis_constant,
            VertexBasedCellPopulation<2>& cell_population
            )
    {
        for (auto cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            // Set cell death threshold
            cell_iter->GetCellData()->SetItem("death threshold", death_threshold);

            // Get SRN model
            DeathClockSrnModel* p_srn_model =
                dynamic_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
            assert(p_srn_model);

            // Set parameters of SRN model
            p_srn_model->SetBaseRate(base_rate);
            p_srn_model->SetNeighboursConstant(neighbours_constant);
            p_srn_model->SetNeighboursG2Constant(neighbours_g2_constant);
            p_srn_model->SetNeighboursApoptosisConstant(neighbours_apoptosis_constant);
        }
    }

public:

    void TestDeathClockBaseRate()
    {
        // Define parameters
        double g1_duration = 10.0;
        double g2_duration = 1e-12;
        double death_threshold = 20.0;
        double base_rate = 2.0;
        double neighbours_constant = 0.0;
        double neighbours_g2_constant = 0.0;
        double neighbours_apoptosis_constant = 0.0;

        // Define birth times
        std::vector<double> birth_times(4, 0.0);

        // Define initial taus
        std::vector<double> initial_taus = { 0.0, 0.0, 0.0, 10.0 };

        // Use the honeycomb vertex mesh generator to create a vertex mesh.
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells =
            create_cells<FixedG1GenerationalCellCycleModel>(
                g1_duration,
                g2_duration,
                birth_times,
                initial_taus);

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Add cell writers
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Define simulation
        OffLatticeSimulation<2> simulator(cell_population);

        // Set death clock parameters
        set_death_clock_parameters(
            death_threshold,
            base_rate,
            neighbours_constant,
            neighbours_g2_constant,
            neighbours_apoptosis_constant,
            cell_population);

        // Configure simulator
        simulator.SetOutputDirectory("TestDeathClockBaseRate");
        simulator.SetSamplingTimestepMultiple(50);

        // Add force
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Add modifiers
        MAKE_PTR(DeathClockModifier<2>, p_death_clock_modifier);
        simulator.AddSimulationModifier(p_death_clock_modifier);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add cell killers
        MAKE_PTR_ARGS(DeathClockCellKiller<2>, p_death_clock_cell_killer,
                (&cell_population));
        simulator.AddCellKiller(p_death_clock_cell_killer);

        // Run simulation for 5 hours
        simulator.SetEndTime(5.0);
        simulator.Solve();

        // Check that taus are as expected
        auto cell_iter = cell_population.Begin();
        auto p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 10.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 10.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 10.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 20.0, 1e-4);

        // Run simulation for another half hour
        simulator.SetEndTime(5.5);
        simulator.Solve();

        // Check that remaining number of cells is 4, with 1 apoptotic cell
        unsigned num_cells = 0;
        unsigned num_apoptotic_cells = 0;
        for (auto cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            num_cells++;
            if (cell_iter->HasApoptosisBegun())
            {
                num_apoptotic_cells++;
            }
        }
        TS_ASSERT_EQUALS(num_cells, 4);
        TS_ASSERT_EQUALS(num_apoptotic_cells, 1);
    }

    void TestDeathClockNeighbours()
    {
        // Define parameters
        double g1_duration = 7.0;
        double g2_duration = 1e-12;
        double death_threshold = 61.0;
        double base_rate = 0.0;
        double neighbours_constant = 2.0;
        double neighbours_g2_constant = 0.0;
        double neighbours_apoptosis_constant = 0.0;

        // Define birth times
        std::vector<double> birth_times(9, 0.0);

        // Define initial taus
        std::vector<double> initial_taus(9, 0.0);

        // Use the honeycomb vertex mesh generator to create a vertex mesh.
        HoneycombVertexMeshGenerator generator(3, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells =
            create_cells<FixedG1GenerationalCellCycleModel>(
                g1_duration,
                g2_duration,
                birth_times,
                initial_taus);

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Add cell writers
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Define simulation
        OffLatticeSimulation<2> simulator(cell_population);

        // Set death clock parameters
        set_death_clock_parameters(
            death_threshold,
            base_rate,
            neighbours_constant,
            neighbours_g2_constant,
            neighbours_apoptosis_constant,
            cell_population);

        // Configure simulator
        simulator.SetOutputDirectory("TestDeathClockNeighbours");
        simulator.SetSamplingTimestepMultiple(50);

        // Add force
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Add modifiers
        MAKE_PTR(DeathClockModifier<2>, p_death_clock_modifier);
        simulator.AddSimulationModifier(p_death_clock_modifier);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add cell killers
        MAKE_PTR_ARGS(DeathClockCellKiller<2>, p_death_clock_cell_killer,
                (&cell_population));
        simulator.AddCellKiller(p_death_clock_cell_killer);

        // Run simulation for 5 hours
        simulator.SetEndTime(5.0);
        simulator.Solve();

        // Check that taus are as expected
        auto cell_iter = cell_population.Begin();
        auto p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 20.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 40.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 30.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 50.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 60.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 30.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 20.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 40.0, 1e-4);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 30.0, 1e-4);

        // Run simulation for another half hour
        simulator.SetEndTime(5.5);
        simulator.Solve();

        // Check that remaining number of cells is 9, with one apoptotic cell
        unsigned num_cells = 0;
        unsigned num_apoptotic_cells = 0;
        for (auto cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            num_cells++;
            if (cell_iter->HasApoptosisBegun())
            {
                num_apoptotic_cells++;
            }
        }
        TS_ASSERT_EQUALS(num_cells, 9);
        TS_ASSERT_EQUALS(num_apoptotic_cells, 1);
    }

    void TestDeathClockNeighboursG2()
    {
        // Define parameters
        double g1_duration = 10.0;
        double g2_duration = 10.0;
        double death_threshold = 30.0;
        double base_rate = 0.0;
        double neighbours_constant = 0.0;
        double neighbours_g2_constant = 2.0;
        double neighbours_apoptosis_constant = 0.0;

        // Define birth times
        std::vector<double> birth_times = { -10.0, 0.0, -10.0, -10.0 };

        // Define initial taus
        std::vector<double> initial_taus(4, 0.0);

        // Use the honeycomb vertex mesh generator to create a vertex mesh.
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells =
            create_cells<FixedG1GenerationalCellCycleModel>(
                g1_duration,
                g2_duration,
                birth_times,
                initial_taus);

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Add cell writers
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Define simulation
        OffLatticeSimulation<2> simulator(cell_population);

        // Set death clock parameters
        set_death_clock_parameters(
            death_threshold,
            base_rate,
            neighbours_constant,
            neighbours_g2_constant,
            neighbours_apoptosis_constant,
            cell_population);

        // Configure simulator
        simulator.SetOutputDirectory("TestDeathClockNeighboursG2");
        simulator.SetSamplingTimestepMultiple(50);

        // Add force
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Add modifiers
        MAKE_PTR(DeathClockModifier<2>, p_death_clock_modifier);
        simulator.AddSimulationModifier(p_death_clock_modifier);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add cell killers
        MAKE_PTR_ARGS(DeathClockCellKiller<2>, p_death_clock_cell_killer,
                (&cell_population));
        simulator.AddCellKiller(p_death_clock_cell_killer);

        // Run simulation for 5 hours
        simulator.SetEndTime(5.0);
        simulator.Solve();

        // Check that taus are as expected
        auto cell_iter = cell_population.Begin();
        auto p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 10.0, 0.1);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 30.0, 0.1);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 20.0, 0.1);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 10.0, 0.1);

        // Run simulation for another half hour
        simulator.SetEndTime(5.5);
        simulator.Solve();

        // Check that remaining number of cells is 4, with one apoptotic cell
        unsigned num_cells = 0;
        unsigned num_apoptotic_cells = 0;
        for (auto cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            num_cells++;
            if (cell_iter->HasApoptosisBegun())
            {
                num_apoptotic_cells++;
            }
        }
        TS_ASSERT_EQUALS(num_cells, 4);
        TS_ASSERT_EQUALS(num_apoptotic_cells, 1);
    }

    void TestDeathClockNeighboursApoptosis()
    {
        // Define parameters
        double g1_duration = 20.0;
        double g2_duration = 1e-12;
        double death_threshold = 10.0;
        double base_rate = 2.0;
        double neighbours_constant = 0.0;
        double neighbours_g2_constant = 0.0;
        double neighbours_apoptosis_constant = -1.0;

        // Define birth times
        std::vector<double> birth_times(4, 0.0);

        // Define initial taus
        std::vector<double> initial_taus = { 0.0, 10.0, 0.0, 0.0 };

        // Use the honeycomb vertex mesh generator to create a vertex mesh.
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells =
            create_cells<FixedG1GenerationalCellCycleModel>(
                g1_duration,
                g2_duration,
                birth_times,
                initial_taus);

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Add cell writers
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        cell_population.AddCellWriter<CellAgesWriter>();

        // Define simulation
        OffLatticeSimulation<2> simulator(cell_population);

        // Set death clock parameters
        set_death_clock_parameters(
            death_threshold,
            base_rate,
            neighbours_constant,
            neighbours_g2_constant,
            neighbours_apoptosis_constant,
            cell_population);

        // Configure simulator
        simulator.SetOutputDirectory("TestDeathClockNeighboursApoptosis");
        simulator.SetSamplingTimestepMultiple(50);

        // Add force
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // Add modifiers
        MAKE_PTR(DeathClockModifier<2>, p_death_clock_modifier);
        simulator.AddSimulationModifier(p_death_clock_modifier);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add cell killers
        MAKE_PTR_ARGS(DeathClockCellKiller<2>, p_death_clock_cell_killer,
                (&cell_population));
        simulator.AddCellKiller(p_death_clock_cell_killer);

        // Run simulation for half an hour
        simulator.SetEndTime(0.5);
        simulator.Solve();

        // Check that taus are as expected
        auto cell_iter = cell_population.Begin();
        auto p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-2);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 10.0, 1e-2);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-2);

        ++cell_iter;
        p_srn_model = static_cast<DeathClockSrnModel*>(cell_iter->GetSrnModel());
        TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-2);

        // Check that remaining number of cells is 4, with one apoptotic cell
        unsigned num_cells = 0;
        unsigned num_apoptotic_cells = 0;
        for (auto cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            num_cells++;
            if (cell_iter->HasApoptosisBegun())
            {
                num_apoptotic_cells++;
            }
        }
        TS_ASSERT_EQUALS(num_cells, 4);
        TS_ASSERT_EQUALS(num_apoptotic_cells, 1);
    }
};

#endif // TESTDEATHCLOCK_HPP_
