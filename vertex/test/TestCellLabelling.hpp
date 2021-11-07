#ifndef TESTCELLLABELLING_HPP_
#define TESTCELLLABELLING_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"

#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "AbstractCellCycleModel.hpp"

#include "CellLabelling.hpp"

class TestCellLabelling : public AbstractCellBasedTestSuite
{

    public:

        void TestRandomlyLabelCells()
        {
            // Set number of cells
            int num_cells = 16;
            double b_cells_ratio = 0.25;

            // Initialise cells
            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellProperty>
                p_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

            // Iteratively generate cells
            for (int i = 0; i < num_cells; i++)
            {
                // Make cell cycle model
                ExponentialG1GenerationalCellCycleModel* p_cell_cycle_model =
                    new ExponentialG1GenerationalCellCycleModel;

                // Wild-type cell mutation state
                boost::shared_ptr<AbstractCellProperty>
                    p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

                // Create cell
                CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

                // Set proliferative type
                p_cell->SetCellProliferativeType(
                        CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

                // Append cell
                cells.push_back(p_cell);
            }

            // Make list of cells
            std::list<CellPtr> cells_list(cells.begin(), cells.end());

            // Randomly label half of the cells
            boost::shared_ptr<AbstractCellProperty>
                p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
            double b_cell_target_area = 1.0;
            double b_cell_g1_duration = 10.0;
            double b_cell_g2_duration = 10.0;
            unsigned b_cell_max_generations = 3;

            RandomlyLabelCells(cells_list, p_state, b_cells_ratio,
                    b_cell_target_area, b_cell_g1_duration, b_cell_g2_duration,
                    b_cell_max_generations);

            // Count number of labelled cells
            unsigned num_labelled_cells = 0;
            for (auto cell_iter = cells.begin();
                    cell_iter != cells.end();
                    ++cell_iter)
            {
                if ((*cell_iter)->template HasCellProperty<CellLabel>())
                {
                    num_labelled_cells++;

                    // Check target area
                    TS_ASSERT_EQUALS((*cell_iter)->GetCellData()->GetItem("target area"), b_cell_target_area);

                    // Check cell cycle
                    ExponentialG1GenerationalCellCycleModel* p_cell_cycle_model =
                        static_cast<ExponentialG1GenerationalCellCycleModel*>(
                                (*cell_iter)->GetCellCycleModel());

                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetTransitCellG1Duration(), b_cell_g1_duration);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetG2Duration(), b_cell_g2_duration);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetMDuration(), 1e-12);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetSDuration(), 1e-12);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetMaxTransitGenerations(), b_cell_max_generations);

                }
            }

            // Check that 4 cells are labelled
            TS_ASSERT_EQUALS(num_labelled_cells, 4);
        }

        void TestDeterministicallyLabelCells()
        {
            // Set number of cells
            int num_cells = 16;
            double b_cells_ratio = 0.25;

            // Initialise cells
            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellProperty>
                p_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

            // Iteratively generate cells
            for (int i = 0; i < num_cells; i++)
            {
                // Make cell cycle model
                ExponentialG1GenerationalCellCycleModel* p_cell_cycle_model =
                    new ExponentialG1GenerationalCellCycleModel;

                // Wild-type cell mutation state
                boost::shared_ptr<AbstractCellProperty>
                    p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

                // Create cell
                CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

                // Set proliferative type
                p_cell->SetCellProliferativeType(
                        CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

                // Append cell
                cells.push_back(p_cell);
            }

            // Make list of cells
            std::list<CellPtr> cells_list(cells.begin(), cells.end());

            // Deterministically label half of the cells
            boost::shared_ptr<AbstractCellProperty>
                p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
            double b_cell_target_area = 1.0;
            double b_cell_g1_duration = 10.0;
            double b_cell_g2_duration = 10.0;
            unsigned b_cell_max_generations = 3;

            DeterministicallyLabelCells(cells_list, p_state, b_cells_ratio,
                    b_cell_target_area, b_cell_g1_duration, b_cell_g2_duration,
                    b_cell_max_generations);

            // Count number of labelled cells
            unsigned num_labelled_cells = 0;
            for (auto cell_iter = cells.begin();
                    cell_iter != cells.end();
                    ++cell_iter)
            {
                if ((*cell_iter)->template HasCellProperty<CellLabel>())
                {
                    num_labelled_cells++;

                    // Check target area
                    TS_ASSERT_EQUALS((*cell_iter)->GetCellData()->GetItem("target area"), b_cell_target_area);

                    // Check cell cycle
                    ExponentialG1GenerationalCellCycleModel* p_cell_cycle_model =
                        static_cast<ExponentialG1GenerationalCellCycleModel*>(
                                (*cell_iter)->GetCellCycleModel());

                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetTransitCellG1Duration(), b_cell_g1_duration);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetG2Duration(), b_cell_g2_duration);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetMDuration(), 1e-12);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetSDuration(), 1e-12);
                    TS_ASSERT_EQUALS(p_cell_cycle_model->GetMaxTransitGenerations(), b_cell_max_generations);
                }

                // First 4 cells should be labelled
                if (num_labelled_cells < 4)
                    TS_ASSERT((*cell_iter)->template HasCellProperty<CellLabel>());
            }

            // Check that 4 cells are labelled
            TS_ASSERT_EQUALS(num_labelled_cells, 4);
        }

};
#endif // TESTCELLLABELLING_HPP_
