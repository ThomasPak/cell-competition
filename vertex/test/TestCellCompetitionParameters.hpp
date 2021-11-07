#ifndef TESTCELLCOMPETITIONPARAMETERS_HPP_
#define TESTCELLCOMPETITIONPARAMETERS_HPP_

#include <cxxtest/TestSuite.h>

#include "CellCompetitionParameters.hpp"

class TestCellCompetitionParameters : public CxxTest::TestSuite
{
    public:

        void TestCellCompetitionParametersDefault()
        {
            CellCompetitionParameters param;

            TS_ASSERT_EQUALS(param.m_output_directory, "cell-competition");
            TS_ASSERT_EQUALS(param.m_simulation_id, "0");

            TS_ASSERT_EQUALS(param.m_seed, 0);

            TS_ASSERT_EQUALS(param.m_simulation_time, 100);
            TS_ASSERT_EQUALS(param.m_dt, 0.01);
            TS_ASSERT_EQUALS(param.m_sampling_timestep_multiple, 100);

            TS_ASSERT_EQUALS(param.m_cell_rearrangement_threshold, 0.01);
            TS_ASSERT_EQUALS(param.m_cell_rearrangement_ratio, 1.5);
            TS_ASSERT_EQUALS(param.m_t2_threshold, 0.001);

            TS_ASSERT_EQUALS(param.m_num_cells_across, 4);
            TS_ASSERT_EQUALS(param.m_num_cells_up, 4);
            TS_ASSERT_EQUALS(param.m_random_labelling, true);
            TS_ASSERT_EQUALS(param.m_initial_b_ratio, 0.5);

            TS_ASSERT_EQUALS(param.m_random_movement_parameter, 0.1);

            TS_ASSERT_EQUALS(param.m_cell_ab_line_tension_parameter, 0.12);

            TS_ASSERT_EQUALS(param.m_cell_a_target_area, 1.0);
            TS_ASSERT_EQUALS(param.m_cell_a_elasticity_parameter, 1.0);
            TS_ASSERT_EQUALS(param.m_cell_a_contractility_parameter, 0.04);
            TS_ASSERT_EQUALS(param.m_cell_a_line_tension_parameter, 0.12);
            TS_ASSERT_EQUALS(param.m_cell_a_boundary_tension_parameter, 0.12);
            TS_ASSERT_EQUALS(param.m_cell_a_g1_duration, 30.0);
            TS_ASSERT_EQUALS(param.m_cell_a_g2_duration, 70.0);
            TS_ASSERT_EQUALS(param.m_cell_a_max_generations, 3);

            TS_ASSERT_EQUALS(param.m_cell_b_target_area, 1.0);
            TS_ASSERT_EQUALS(param.m_cell_b_elasticity_parameter, 1.0);
            TS_ASSERT_EQUALS(param.m_cell_b_contractility_parameter, 0.04);
            TS_ASSERT_EQUALS(param.m_cell_b_line_tension_parameter, 0.12);
            TS_ASSERT_EQUALS(param.m_cell_b_boundary_tension_parameter, 0.12);
            TS_ASSERT_EQUALS(param.m_cell_b_g1_duration, 30.0);
            TS_ASSERT_EQUALS(param.m_cell_b_g2_duration, 70.0);
            TS_ASSERT_EQUALS(param.m_cell_b_max_generations, 3);
        }

        void TestCellCompetitionParametersUpdate1()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("output-directory = other-dir");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_output_directory, "other-dir");
        }

        void TestCellCompetitionParametersUpdate2()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("simulation-time = 50");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_simulation_time, 50);
        }

        void TestCellCompetitionParametersUpdate3()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("num-cells-across = 10");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_num_cells_across, 10);
        }

        void TestCellCompetitionParametersUpdate4()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("num-cells-up = 10");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_num_cells_up, 10);
        }

        void TestCellCompetitionParametersUpdate5()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("random-movement-parameter = 0.5");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_random_movement_parameter, 0.5);
        }

        void TestCellCompetitionParametersUpdate6()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-rearrangement-threshold = 0.05");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_rearrangement_threshold, 0.05);
        }

        void TestCellCompetitionParametersUpdate7()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("dt = 0.0025");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_dt, 0.0025);
        }

        void TestCellCompetitionParametersUpdate8()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("sampling-timestep-multiple = 100");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_sampling_timestep_multiple, 100);
        }

        void TestCellCompetitionParametersUpdate9()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-ab-line-tension-parameter = 0.25");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_ab_line_tension_parameter, 0.25);
        }

        void TestCellCompetitionParametersUpdate10()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-target-area = 0.5");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_target_area, 0.5);
        }

        void TestCellCompetitionParametersUpdate11()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-elasticity-parameter = 50.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_elasticity_parameter, 50.0);
        }

        void TestCellCompetitionParametersUpdate12()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-contractility-parameter = 5.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_contractility_parameter, 5.0);
        }

        void TestCellCompetitionParametersUpdate13()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-line-tension-parameter = 0.25");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_line_tension_parameter, 0.25);
        }

        void TestCellCompetitionParametersUpdate14()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-boundary-tension-parameter = 0.5");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_boundary_tension_parameter, 0.5);
        }

        void TestCellCompetitionParametersUpdate15()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-target-area = 0.5");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_target_area, 0.5);
        }

        void TestCellCompetitionParametersUpdate16()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-elasticity-parameter = 50.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_elasticity_parameter, 50.0);
        }

        void TestCellCompetitionParametersUpdate17()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-contractility-parameter = 5.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_contractility_parameter, 5.0);
        }

        void TestCellCompetitionParametersUpdate18()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-line-tension-parameter = 0.25");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_line_tension_parameter, 0.25);
        }

        void TestCellCompetitionParametersUpdate19()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-boundary-tension-parameter = 0.5");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_boundary_tension_parameter, 0.5);
        }

        void TestCellCompetitionParametersUpdate20()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("initial-b-ratio = 0.25");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_initial_b_ratio, 0.25);
        }

        void TestCellCompetitionParametersUpdate21()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-g1-duration = 1.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_g1_duration, 1.0);
        }

        void TestCellCompetitionParametersUpdate22()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-g1-duration = 1.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_g1_duration, 1.0);
        }

        void TestCellCompetitionParametersUpdate23()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-max-generations = 2");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_max_generations, 2);
        }

        void TestCellCompetitionParametersUpdate24()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-max-generations = 2");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_max_generations, 2);
        }

        void TestCellCompetitionParametersUpdate25()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("simulation-id = 2");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_simulation_id, "2");
        }

        void TestCellCompetitionParametersUpdate26()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("seed = 3");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_seed, 3);
        }

        void TestCellCompetitionParametersUpdate27()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-a-g2-duration = 1.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_a_g2_duration, 1.0);
        }

        void TestCellCompetitionParametersUpdate28()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-b-g2-duration = 1.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_b_g2_duration, 1.0);
        }

        void TestCellCompetitionParametersUpdate29()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("random-labelling = false");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_random_labelling, false);
        }

        void TestCellCompetitionParametersUpdate30()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("cell-rearrangement-ratio = 2.0");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_cell_rearrangement_ratio, 2.0);
        }

        void TestCellCompetitionParametersUpdate31()
        {
            CellCompetitionParameters param;

            Dictionary dict = create_dictionary("t2-threshold = 0.01");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_t2_threshold, 0.01);
        }
};

#endif // TESTCELLCOMPETITIONPARAMETERS_HPP_
