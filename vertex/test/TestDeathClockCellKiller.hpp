#ifndef TESTDEATHCLOCKCELLKILLER_HPP_
#define TESTDEATHCLOCKCELLKILLER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "UniformCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "DeathClockCellKiller.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestDeathClockCellKiller : public AbstractCellBasedTestSuite
{
    public:

        void TestDeathClockCellKiller1()
        {
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

            // Create cell killer
            DeathClockCellKiller<2> death_clock_cell_killer(&cell_population);

            TS_ASSERT_EQUALS(death_clock_cell_killer.GetIdentifier(), "DeathClockCellKiller-2");

            // Set tau of cells to incrementing values and set death threshold
            // to 5.0
            double tau = 0.0;
            double death_threshold = 5.0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem("death threshold",
                        death_threshold);
                cell_iter->GetCellData()->SetItem("tau", tau);
                ++tau;
            }

            // Check that some of the vector of cells reach apotosis
            death_clock_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

            std::set<double> old_locations;

            unsigned num_apoptotic_cells = 0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                if (cell_iter->HasApoptosisBegun())
                {
                    ++num_apoptotic_cells;
                }
            }

            TS_ASSERT_EQUALS(num_apoptotic_cells, 4);

            // DeathClockCellKiller does not set a death time when labelling
            // cells as apoptotic, so that all cells are removed by T2 swaps.
            // Hence, we do not test whether DeathClockCellKiller by itself
            // kills any cells.
        }

        // Test that cells who are in G1 phase are killed off
        void TestDeathClockCellKiller2()
        {
            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a transit cell type and
            // FixedG1GenerationalCellCycleModel
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(TransitCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create cell killer
            DeathClockCellKiller<2> death_clock_cell_killer(&cell_population);

            TS_ASSERT_EQUALS(death_clock_cell_killer.GetIdentifier(), "DeathClockCellKiller-2");

            // Set tau of cells to incrementing values and set death threshold
            // to 5.0.  The birth time of cells have been initialised to
            // different values, such that two cells are in G1 phase and will
            // be labelled as apoptotic
            double tau = 6.0;
            double death_threshold = 5.0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem("death threshold",
                        death_threshold);
                cell_iter->GetCellData()->SetItem("tau", tau);

                cell_iter->GetCellCycleModel()->Initialise();
                cell_iter->ReadyToDivide();
            }

            // Check that some of the vector of cells reach apotosis
            death_clock_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

            std::set<double> old_locations;

            unsigned num_apoptotic_cells = 0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                if (cell_iter->HasApoptosisBegun())
                {
                    ++num_apoptotic_cells;
                }
            }

            TS_ASSERT_EQUALS(num_apoptotic_cells, 2);

            // DeathClockCellKiller does not set a death time when labelling
            // cells as apoptotic, so that all cells are removed by T2 swaps.
            // Hence, we do not test whether DeathClockCellKiller by itself
            // kills any cells.
        }

        // Same as TestDeathClockCellKiller2, but with
        // ExponentialG1GenerationalCellCycleModel
        void TestDeathClockCellKiller3()
        {
            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a transit cell type and
            // ExponentialG1GenerationalCellCycleModel
            CellsGenerator<ExponentialG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(TransitCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create cell killer
            DeathClockCellKiller<2> death_clock_cell_killer(&cell_population);

            TS_ASSERT_EQUALS(death_clock_cell_killer.GetIdentifier(), "DeathClockCellKiller-2");

            // Set tau of cells to incrementing values and set death threshold
            // to 5.0
            double tau = 6.0;
            double death_threshold = 5.0;
            unsigned num_cells_g1 = 0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem("death threshold",
                        death_threshold);
                cell_iter->GetCellData()->SetItem("tau", tau);

                ExponentialG1GenerationalCellCycleModel* p_ccm =
                    dynamic_cast<ExponentialG1GenerationalCellCycleModel*>(
                            cell_iter->GetCellCycleModel());
                assert(p_ccm);

                p_ccm->SetRate(1.0 / 10.0);

                cell_iter->GetCellCycleModel()->Initialise();
                cell_iter->ReadyToDivide();

                if (p_ccm->GetCurrentCellCyclePhase() == G_ONE_PHASE)
                {
                    ++num_cells_g1;
                }
            }

            // Check that some of the vector of cells reach apotosis
            death_clock_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

            std::set<double> old_locations;

            unsigned num_apoptotic_cells = 0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                if (cell_iter->HasApoptosisBegun())
                {
                    ++num_apoptotic_cells;
                }
            }

            TS_ASSERT_EQUALS(num_apoptotic_cells, num_cells_g1);

            // DeathClockCellKiller does not set a death time when labelling
            // cells as apoptotic, so that all cells are removed by T2 swaps.
            // Hence, we do not test whether DeathClockCellKiller by itself
            // kills any cells.
        }
};
#endif // TESTDEATHCLOCKCELLKILLER_HPP_
