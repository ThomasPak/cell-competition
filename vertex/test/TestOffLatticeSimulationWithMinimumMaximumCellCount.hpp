#ifndef TESTOFFLATTICESIMULATIONWITHCELLCOUNTLIMIT_HPP
#define TESTOFFLATTICESIMULATIONWITHCELLCOUNTLIMIT_HPP

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <cmath>

#include <unordered_map>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CustomUniformG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "OffLatticeSimulationWithMinimumMaximumCellCount.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "LogFile.hpp"
#include "TargetedCellKiller.hpp"
#include "CellAncestor.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestOffLatticeSimulation : public AbstractCellBasedWithTimingsTestSuite
{
public:
    void TestOffLatticeSimulationWithMaximumCellCount()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        unsigned max_cell_count = 7;

        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Define mutation state and proliferative type
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Create cells
        std::vector<CellPtr> cells;
        std::vector<double> birth_times = { -10.0, -13.2, -50.0, -50.0 };
        for (const auto& birth_time : birth_times)
        {
            // Get individual parameters
            double g1_duration = 50.0;
            double g2_duration = 50.0;
            double uniform_g1_range = 0.0;
            unsigned max_transit_generation = static_cast<unsigned>(-1);

            // Create cell cycle model
            CustomUniformG1GenerationalCellCycleModel* p_cell_cycle_model =
                new CustomUniformG1GenerationalCellCycleModel();

            // Set cell cycle phase durations
            p_cell_cycle_model->SetMDuration(1e-12);
            p_cell_cycle_model->SetSDuration(1e-12);

            p_cell_cycle_model->SetTransitCellG1Duration(g1_duration);
            p_cell_cycle_model->SetG2Duration(g2_duration);

            p_cell_cycle_model->SetRange(uniform_g1_range);

            // Set max transit generations
            p_cell_cycle_model->SetMaxTransitGenerations(max_transit_generation);

            // Create cell
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            // Set cell proliferative type
            p_cell->SetCellProliferativeType(p_transit_type);

            // Set cell birth time
            p_cell->SetBirthTime(birth_time);

            // Push cell
            cells.push_back(p_cell);
        }

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Define simulation
        OffLatticeSimulationWithMinimumMaximumCellCount simulator(cell_population, 0, max_cell_count);
        simulator.SetOutputDirectory("TestCellPopulationSimWithMaximumCellCount");

        TS_ASSERT_EQUALS(simulator.GetMaximumCellCount(), max_cell_count);

        // Set the end time to 100h, the stopping event occurs at t = 86.8
        // hours however
        simulator.SetEndTime(100.0);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);


        // Run cell-based simulation
        simulator.Solve();

        double time = SimulationTime::Instance()->GetTime();
        TS_ASSERT_DELTA(time, 87.0, 1e-1);

        // The number of cells should be greater than or equal to the maximum
        // cell count
        unsigned cell_count = 0;
        for (auto cell_ptr = simulator.rGetCellPopulation().Begin();
                cell_ptr != simulator.rGetCellPopulation().End();
                ++cell_ptr)
        {
            cell_count++;
        }

        TS_ASSERT_LESS_THAN_EQUALS(max_cell_count, cell_count);
    }

    void TestOffLatticeSimulationWithMinimumCellCount()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        unsigned min_cell_count = 1;

        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Define mutation state and proliferative type
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Create cells
        std::vector<CellPtr> cells;
        std::vector<double> birth_times = { -10.0, -13.2, -50.0, -50.0 };
        for (const auto& birth_time : birth_times)
        {
            // Get individual parameters
            double g1_duration = 50.0;
            double g2_duration = 50.0;
            double uniform_g1_range = 0.0;
            unsigned max_transit_generation = static_cast<unsigned>(-1);

            // Create cell cycle model
            CustomUniformG1GenerationalCellCycleModel* p_cell_cycle_model =
                new CustomUniformG1GenerationalCellCycleModel();

            // Set cell cycle phase durations
            p_cell_cycle_model->SetMDuration(1e-12);
            p_cell_cycle_model->SetSDuration(1e-12);

            p_cell_cycle_model->SetTransitCellG1Duration(g1_duration);
            p_cell_cycle_model->SetG2Duration(g2_duration);

            p_cell_cycle_model->SetRange(uniform_g1_range);

            // Set max transit generations
            p_cell_cycle_model->SetMaxTransitGenerations(max_transit_generation);

            // Create cell
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            // Set cell proliferative type
            p_cell->SetCellProliferativeType(p_transit_type);

            // Set cell birth time
            p_cell->SetBirthTime(birth_time);

            // Push cell
            cells.push_back(p_cell);
        }

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Define simulation
        OffLatticeSimulationWithMinimumMaximumCellCount
            simulator(cell_population, min_cell_count);
        simulator.SetOutputDirectory("TestCellPopulationSimWithMinimumCellCount");

        TS_ASSERT_EQUALS(simulator.GetMinimumCellCount(), min_cell_count);

        // Set the end time to 100h, the stopping event occurs at t = 0.01
        // hours however
        simulator.SetEndTime(100.0);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);

        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer_1, (&cell_population, 3u));
        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer_2, (&cell_population, 2u));
        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer_3, (&cell_population, 1u));

        simulator.AddCellKiller(p_killer_1);
        simulator.AddCellKiller(p_killer_2);
        simulator.AddCellKiller(p_killer_3);

        // Run cell-based simulation
        simulator.Solve();

        double time = SimulationTime::Instance()->GetTime();
        TS_ASSERT_DELTA(time, 1, 1e-1);

        // The number of cells should be less or equal to the cell count limit
        unsigned cell_count = 0;
        for (auto cell_ptr = simulator.rGetCellPopulation().Begin();
                cell_ptr != simulator.rGetCellPopulation().End();
                ++cell_ptr)
        {
            cell_count++;
        }

        TS_ASSERT_LESS_THAN_EQUALS(cell_count, min_cell_count);
    }

    void TestOffLatticeSimulationWithMaximumCellCountForAncestor()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        std::unordered_map<unsigned, unsigned> min_cell_count_for_ancestor;
        std::unordered_map<unsigned, unsigned> max_cell_count_for_ancestor;
        max_cell_count_for_ancestor[0] = 3;

        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Define mutation state and proliferative type
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Create cells
        std::vector<CellPtr> cells;
        std::vector<double> birth_times = { -10.0, -13.2, -50.0, -50.0 };
        for (const auto& birth_time : birth_times)
        {
            // Get individual parameters
            double g1_duration = 50.0;
            double g2_duration = 50.0;
            double uniform_g1_range = 0.0;
            unsigned max_transit_generation = static_cast<unsigned>(-1);

            // Create cell cycle model
            CustomUniformG1GenerationalCellCycleModel* p_cell_cycle_model =
                new CustomUniformG1GenerationalCellCycleModel();

            // Set cell cycle phase durations
            p_cell_cycle_model->SetMDuration(1e-12);
            p_cell_cycle_model->SetSDuration(1e-12);

            p_cell_cycle_model->SetTransitCellG1Duration(g1_duration);
            p_cell_cycle_model->SetG2Duration(g2_duration);

            p_cell_cycle_model->SetRange(uniform_g1_range);

            // Set max transit generations
            p_cell_cycle_model->SetMaxTransitGenerations(max_transit_generation);

            // Create cell
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            // Set cell proliferative type
            p_cell->SetCellProliferativeType(p_transit_type);

            // Set cell birth time
            p_cell->SetBirthTime(birth_time);

            // Push cell
            cells.push_back(p_cell);
        }

        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor_0, (0));
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor_1, (1));

        cells[0]->SetAncestor(p_cell_ancestor_0);
        cells[1]->SetAncestor(p_cell_ancestor_0);
        cells[2]->SetAncestor(p_cell_ancestor_1);
        cells[3]->SetAncestor(p_cell_ancestor_1);

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Define simulation
        OffLatticeSimulationWithMinimumMaximumCellCount
            simulator(cell_population, 0, static_cast<unsigned>(-1),
                    min_cell_count_for_ancestor, max_cell_count_for_ancestor);

        simulator.SetOutputDirectory("TestCellPopulationSimWithMaximumCellCountForAncestor");

        TS_ASSERT_EQUALS(simulator.GetMaximumCellCountForAncestor(0),
                max_cell_count_for_ancestor[0]);

        // Set the end time to 100h, the stopping event occurs at t = 86.8
        // hours however
        simulator.SetEndTime(100.0);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);

        // Run cell-based simulation
        simulator.Solve();

        double time = SimulationTime::Instance()->GetTime();
        TS_ASSERT_DELTA(time, 87.0, 1e-1);

        // The number of cells of ancestor 0 should be greater than or equal to
        // the maximum cell count
        unsigned cell_count_for_ancestor_0 = 0;
        for (auto cell_ptr = simulator.rGetCellPopulation().Begin();
                cell_ptr != simulator.rGetCellPopulation().End();
                ++cell_ptr)
        {
            if (cell_ptr->GetAncestor() == 0)
            {
                cell_count_for_ancestor_0++;
            }
        }

        TS_ASSERT_LESS_THAN_EQUALS(max_cell_count_for_ancestor[0], cell_count_for_ancestor_0);
    }

    void TestOffLatticeSimulationWithMinimumCellCountForAncestor()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        std::unordered_map<unsigned, unsigned> min_cell_count_for_ancestor;
        std::unordered_map<unsigned, unsigned> max_cell_count_for_ancestor;
        min_cell_count_for_ancestor[1] = 0;

        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Define mutation state and proliferative type
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Create cells
        std::vector<CellPtr> cells;
        std::vector<double> birth_times = { -10.0, -13.2, -50.0, -50.0 };
        for (const auto& birth_time : birth_times)
        {
            // Get individual parameters
            double g1_duration = 50.0;
            double g2_duration = 50.0;
            double uniform_g1_range = 0.0;
            unsigned max_transit_generation = static_cast<unsigned>(-1);

            // Create cell cycle model
            CustomUniformG1GenerationalCellCycleModel* p_cell_cycle_model =
                new CustomUniformG1GenerationalCellCycleModel();

            // Set cell cycle phase durations
            p_cell_cycle_model->SetMDuration(1e-12);
            p_cell_cycle_model->SetSDuration(1e-12);

            p_cell_cycle_model->SetTransitCellG1Duration(g1_duration);
            p_cell_cycle_model->SetG2Duration(g2_duration);

            p_cell_cycle_model->SetRange(uniform_g1_range);

            // Set max transit generations
            p_cell_cycle_model->SetMaxTransitGenerations(max_transit_generation);

            // Create cell
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            // Set cell proliferative type
            p_cell->SetCellProliferativeType(p_transit_type);

            // Set cell birth time
            p_cell->SetBirthTime(birth_time);

            // Push cell
            cells.push_back(p_cell);
        }

        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor_0, (0));
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor_1, (1));

        cells[0]->SetAncestor(p_cell_ancestor_0);
        cells[1]->SetAncestor(p_cell_ancestor_0);
        cells[2]->SetAncestor(p_cell_ancestor_1);
        cells[3]->SetAncestor(p_cell_ancestor_1);

        // Define cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Define simulation
        OffLatticeSimulationWithMinimumMaximumCellCount
            simulator(cell_population, 0, static_cast<unsigned>(-1),
                    min_cell_count_for_ancestor, max_cell_count_for_ancestor);
        simulator.SetOutputDirectory("TestCellPopulationSimWithMinimumCellCountForAncestor");

        TS_ASSERT_EQUALS(simulator.GetMinimumCellCountForAncestor(1),
                min_cell_count_for_ancestor[1]);

        // Set the end time to 100h, the stopping event occurs at t = 0.01
        // hours however
        simulator.SetEndTime(100.0);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);

        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer_1, (&cell_population, 3u));
        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer_2, (&cell_population, 2u));
        MAKE_PTR_ARGS(TargetedCellKiller<2>, p_killer_3, (&cell_population, 1u));

        simulator.AddCellKiller(p_killer_1);
        simulator.AddCellKiller(p_killer_2);
        simulator.AddCellKiller(p_killer_3);

        // Run cell-based simulation
        simulator.Solve();

        double time = SimulationTime::Instance()->GetTime();
        TS_ASSERT_DELTA(time, 1, 1e-1);

        // The number of cells of ancestor 1 should be less or equal to the
        // cell count limit
        unsigned cell_count_for_ancestor_1 = 0;
        for (auto cell_ptr = simulator.rGetCellPopulation().Begin();
                cell_ptr != simulator.rGetCellPopulation().End();
                ++cell_ptr)
        {
            if (cell_ptr->GetAncestor() == 1)
            {
                cell_count_for_ancestor_1++;
            }
        }

        TS_ASSERT_LESS_THAN_EQUALS(cell_count_for_ancestor_1, min_cell_count_for_ancestor[1]);
    }
};

#endif // TESTOFFLATTICESIMULATIONWITHCELLCOUNTLIMIT_HPP
