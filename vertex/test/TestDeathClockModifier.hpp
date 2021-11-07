#ifndef TESTDEATHCLOCKMODIFIER_HPP_
#define TESTDEATHCLOCKMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"

#include "DeathClockModifier.hpp"
#include "DeathClockCellKiller.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "UniformCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

static constexpr double EDGE_LENGTH = pow(3.0, -0.5);

class TestDeathClockModifier : public AbstractCellBasedTestSuite
{
    public:

        void TestNumberOfNeighbours()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestDeathClockModifier1");
            simulator.SetEndTime(0.01);

            // Add Death Clock modifier
            MAKE_PTR(DeathClockModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Check that number of neighbours is correct
            auto cell_iter = cell_population.Begin();

            // First cell should have two neighbours
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 2.0);

            // Second cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 4.0);

            // Third cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);

            // Fourth cell should have five neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 5.0);

            // Fifth cell should have six neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 6.0);

            // Sixth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);

            // Seventh cell should have two neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 2.0);

            // Eighth cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 4.0);

            // Ninth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);
        }

        void TestNumberOfNeighboursG2()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(TransitCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestDeathClockModifier2");
            simulator.SetEndTime(2.5);

            // Add Death Clock modifier
            MAKE_PTR(DeathClockModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Check that number of neighbours is correct
            auto cell_iter = cell_population.Begin();

            // First cell should have two neighbours
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 2.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 0.0);

            // Second cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 4.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 0.0);

            // Third cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 0.0);

            // Fourth cell should have five neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 5.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 2.0);

            // Fifth cell should have six neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 6.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 2.0);

            // Sixth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 1.0);

            // Seventh cell should have two neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 2.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 1.0);

            // Eighth cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 4.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 2.0);

            // Ninth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours g2"), 1.0);
        }

        void TestNumberOfNeighboursApoptosis()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(TransitCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestDeathClockModifier3");
            simulator.SetEndTime(0.1);

            // Add Death Clock modifier
            MAKE_PTR(DeathClockModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

			// Add DeathClockCellKiller
            double death_threshold = 10.0;
            unsigned idx = 0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem("death threshold", death_threshold);

                if (idx == 0 || idx == 1 || idx == 6)
                {
                    cell_iter->SetBirthTime(-1.0);
                    cell_iter->GetCellData()->SetItem("tau", 11.0);
                }
                else
                {
                    cell_iter->GetCellData()->SetItem("tau", 0.0);
                }
                idx++;
            }

			MAKE_PTR_ARGS(DeathClockCellKiller<2>, p_death_clock_cell_killer,
					(&cell_population));
			simulator.AddCellKiller(p_death_clock_cell_killer);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Check that number of neighbours is correct
            auto cell_iter = cell_population.Begin();

            // First cell should have two neighbours
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 2.0);
            TS_ASSERT_EQUALS(cell_iter->HasApoptosisBegun(), true);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 1.0);

            // Second cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 4.0);
            TS_ASSERT_EQUALS(cell_iter->HasApoptosisBegun(), true);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 1.0);

            // Third cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 1.0);

            // Fourth cell should have five neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 5.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 3.0);

            // Fifth cell should have six neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 6.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 1.0);

            // Sixth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 0.0);

            // Seventh cell should have two neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 2.0);
            TS_ASSERT_EQUALS(cell_iter->HasApoptosisBegun(), true);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 0.0);

            // Eighth cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 4.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 1.0);

            // Ninth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours"), 3.0);
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem("number of neighbours apoptosis"), 0.0);
        }

        void TestTotalEdgeLength()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestDeathClockModifier4");
            simulator.SetEndTime(0.001);

            // Add Death Clock modifier
            MAKE_PTR(DeathClockModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Check that number of neighbours is correct
            auto cell_iter = cell_population.Begin();

            // First cell should have two neighbours
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Second cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Third cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Fourth cell should have five neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Fifth cell should have six neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Sixth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Seventh cell should have two neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Eighth cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Ninth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
        }

        void TestSharedEdgeLength()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestDeathClockModifier5");
            simulator.SetEndTime(0.001);

            // Add Death Clock modifier
            MAKE_PTR(DeathClockModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Check that number of neighbours is correct
            auto cell_iter = cell_population.Begin();

            // First cell should have two neighbours
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 2.0 * EDGE_LENGTH, 1e-5);

            // Second cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 4.0 * EDGE_LENGTH, 1e-5);

            // Third cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 3.0 * EDGE_LENGTH, 1e-5);

            // Fourth cell should have five neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 5.0 * EDGE_LENGTH, 1e-5);

            // Fifth cell should have six neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 6.0 * EDGE_LENGTH, 1e-5);

            // Sixth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 3.0 * EDGE_LENGTH, 1e-5);

            // Seventh cell should have two neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 2.0 * EDGE_LENGTH, 1e-5);

            // Eighth cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 4.0 * EDGE_LENGTH, 1e-5);

            // Ninth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length"), 3.0 * EDGE_LENGTH, 1e-5);
        }

        void TestSharedEdgeLengthG2()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(TransitCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestDeathClockModifier6");
            simulator.SetEndTime(2.5);

            // Add Death Clock modifier
            MAKE_PTR(DeathClockModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Check that number of neighbours is correct
            auto cell_iter = cell_population.Begin();

            // First cell should have two neighbours
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 0.0 * EDGE_LENGTH, 1e-5);

            // Second cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 0.0 * EDGE_LENGTH, 1e-5);

            // Third cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 0.0 * EDGE_LENGTH, 1e-5);

            // Fourth cell should have five neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 2.0 * EDGE_LENGTH, 1e-5);

            // Fifth cell should have six neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 2.0 * EDGE_LENGTH, 1e-5);

            // Sixth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 1.0 * EDGE_LENGTH, 1e-5);

            // Seventh cell should have two neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 1.0 * EDGE_LENGTH, 1e-5);

            // Eighth cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 2.0 * EDGE_LENGTH, 1e-5);

            // Ninth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length g2"), 1.0 * EDGE_LENGTH, 1e-5);
        }

        void TestSharedEdgeLengthApoptosis()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(TransitCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("TestDeathClockModifier7");
            simulator.SetEndTime(0.1);

            // Add Death Clock modifier
            MAKE_PTR(DeathClockModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

			// Add DeathClockCellKiller
            double death_threshold = 10.0;
            unsigned idx = 0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem("death threshold", death_threshold);

                if (idx == 0 || idx == 1 || idx == 6)
                {
                    cell_iter->SetBirthTime(-1.0);
                    cell_iter->GetCellData()->SetItem("tau", 11.0);
                }
                else
                {
                    cell_iter->GetCellData()->SetItem("tau", 0.0);
                }
                idx++;
            }

			MAKE_PTR_ARGS(DeathClockCellKiller<2>, p_death_clock_cell_killer,
					(&cell_population));
			simulator.AddCellKiller(p_death_clock_cell_killer);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Check that number of neighbours is correct
            auto cell_iter = cell_population.Begin();

            // First cell should have two neighbours
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->HasApoptosisBegun(), true, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 1.0 * EDGE_LENGTH, 1e-5);

            // Second cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->HasApoptosisBegun(), true, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 1.0 * EDGE_LENGTH, 1e-5);

            // Third cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 1.0 * EDGE_LENGTH, 1e-5);

            // Fourth cell should have five neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 3.0 * EDGE_LENGTH, 1e-5);

            // Fifth cell should have six neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 1.0 * EDGE_LENGTH, 1e-5);

            // Sixth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 0.0 * EDGE_LENGTH, 1e-5);

            // Seventh cell should have two neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->HasApoptosisBegun(), true, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 0.0 * EDGE_LENGTH, 1e-5);

            // Eighth cell should have four neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 1.0 * EDGE_LENGTH, 1e-5);

            // Ninth cell should have three neighbours
            ++cell_iter;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("total edge length"), 6.0 * EDGE_LENGTH, 1e-5);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("shared edge length apoptosis"), 0.0 * EDGE_LENGTH, 1e-5);
        }
};

#endif // TESTDEATHCLOCKMODIFIER_HPP_
