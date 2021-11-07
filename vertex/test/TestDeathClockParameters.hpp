#ifndef TESTDEATHCLOCKPARAMETERS_HPP
#define TESTDEATHCLOCKPARAMETERS_HPP

#include <cxxtest/TestSuite.h>

#include "DeathClockParameters.hpp"

# define TS_ASSERT_EQUALS_VECTORS(vec1, vec2) { \
    TS_ASSERT_EQUALS(vec1.size(), vec2.size()); \
    for (unsigned i = 0; i < vec1.size(); i++) \
    { \
        TS_ASSERT_EQUALS(vec1[i], vec2[i]); \
    } \
}

# define TS_ASSERT_LESS_THAN_EQUALS_VECTORS(vec1, vec2) { \
    TS_ASSERT_EQUALS(vec1.size(), vec2.size()); \
    for (unsigned i = 0; i < vec1.size(); i++) \
    { \
        TS_ASSERT_LESS_THAN_EQUALS(vec1[i], vec2[i]); \
    } \
}

// Hack: include `catch` keyword to enable exception handling
class TestDeathClockParameters : public CxxTest::TestSuite
{
    public:

        void TestDeathClockParametersReadDoubleVector()
        {
            // Test ReadDoubleVector
            std::string input("2.3,4.3,3.3,2.2,2.3,2.4,2.3,-1.3,5.3,2.0,2.5,2.1,2.1,2.6,5.5,3.3");
            std::vector<double> test_output = {2.3,4.3,3.3,2.2,2.3,2.4,2.3,-1.3,5.3,2.0,2.5,2.1,2.1,2.6,5.5,3.3};
            std::vector<double> output = DeathClockParameters::ReadDoubleVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Nonsense should be rejected
            input = "afds re2, 434 434889324";
            TS_ASSERT_THROWS_ANYTHING(DeathClockParameters::ReadDoubleVector(input));
        }

        void TestDeathClockParametersReadUnsignedVector()
        {
            // Test ReadUnsignedVector
            std::string input("2,4,3,2,2,2,2,3,5,2,2,2,2,2,5,3");
            std::vector<unsigned> test_output = {2,4,3,2,2,2,2,3,5,2,2,2,2,2,5,3};
            std::vector<unsigned> output = DeathClockParameters::ReadUnsignedVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Nonsense should be rejected
            input = "afds re2, 434 434889324";
            TS_ASSERT_THROWS_ANYTHING(DeathClockParameters::ReadUnsignedVector(input));
        }

        void TestDeathClockParametersReadUnsignedToUnsignedMap()
        {
            // Test ReadUnsignedToUnsignedMap
            std::string input("{0:1, 5 : 3, 3:  7}");
            std::map<unsigned, unsigned> test_output;
            test_output[0] = 1;
            test_output[5] = 3;
            test_output[3] = 7;
            auto output =
                DeathClockParameters::ReadUnsignedToUnsignedMap(input);

            for (const auto& item : output)
            {
                TS_ASSERT_EQUALS(output.at(item.first), test_output.at(item.first));
            }

            // Nonsense should be rejected
            input = "afds re2, 434 434889324";
            TS_ASSERT_THROWS_ANYTHING(DeathClockParameters::ReadUnsignedToUnsignedMap(input));
        }

        void TestDeathClockParametersPrintUnsignedToUnsignedMap()
        {
            // Test PrintUnsignedToUnsignedMap
            std::map<unsigned, unsigned> input;
            input[0] = 1;
            input[5] = 3;
            input[3] = 7;
            auto output =
                DeathClockParameters::PrintUnsignedToUnsignedMap(input);

            TS_ASSERT_EQUALS(output, "{0:1, 3:7, 5:3}");
            TS_ASSERT_DIFFERS(output, "{0:1, 5:3, 3:7}");
        }

        void TestDeathClockParametersProcessDoubleVector()
        {
            DeathClockParameters param;

            // Nonsense should be rejected
            std::string input = "afds re2, 434 434889324";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessDoubleVector(input));

            // Test ProcessDoubleVector with vector size == m_initial_num_cells
            input = "2.3,4.3,3.3,2.2,2.3,2.4,2.3,-1.3,5.3,2.0,2.5,2.1,2.1,2.6,5.5,3.3";
            std::vector<double> test_output = {2.3,4.3,3.3,2.2,2.3,2.4,2.3,-1.3,5.3,2.0,2.5,2.1,2.1,2.6,5.5,3.3};
            std::vector<double> output = param.ProcessDoubleVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Test ProcessDoubleVector with vector size == 1
            input = "2.3";
            test_output = std::vector<double>(param.m_initial_num_cells, 2.3);
            output = param.ProcessDoubleVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Other vector sizes should raise an error
            input = "2.3,5";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessDoubleVector(input));
            input = "";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessDoubleVector(input));

            // Try the same with m_initial_num_cells = 2
            // Note that this is no longer a consistent set of parameters
            param.m_initial_num_cells = 2;
            param.m_num_cells_across = 2;
            param.m_num_cells_up = 1;

            // Test ProcessDoubleVector with vector size == m_initial_num_cells
            input = "2.3,4.3";
            test_output = {2.3,4.3};
            output = param.ProcessDoubleVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Test ProcessDoubleVector with vector size == 1
            input = "2.3";
            test_output = std::vector<double>(param.m_initial_num_cells, 2.3);
            output = param.ProcessDoubleVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Other vector sizes should raise an error
            input = "2.3,5,2";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessDoubleVector(input));
            input = "";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessDoubleVector(input));
        }

        void TestDeathClockParametersProcessUnsignedVector()
        {
            DeathClockParameters param;

            // Nonsense should be rejected
            std::string input = "afds re2, 434 434889324";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessUnsignedVector(input));

            // Test ProcessUnsignedVector with vector size == m_initial_num_cells
            input = "2,4,3,2,2,2,2,1,5,2,2,2,2,2,5,3";
            std::vector<unsigned> test_output = {2,4,3,2,2,2,2,1,5,2,2,2,2,2,5,3};
            std::vector<unsigned> output = param.ProcessUnsignedVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Test ProcessUnsignedVector with vector size == 1
            input = "2";
            test_output = std::vector<unsigned>(param.m_initial_num_cells, 2);
            output = param.ProcessUnsignedVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Other vector sizes should raise an error
            input = "2,5";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessUnsignedVector(input));
            input = "";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessUnsignedVector(input));

            // Try the same with m_initial_num_cells = 2
            // Note that this is no longer a consistent set of parameters
            param.m_initial_num_cells = 2;
            param.m_num_cells_across = 2;
            param.m_num_cells_up = 1;

            // Test ProcessUnsignedVector with vector size == m_initial_num_cells
            input = "2,4";
            test_output = {2,4};
            output = param.ProcessUnsignedVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Test ProcessUnsignedVector with vector size == 1
            input = "2";
            test_output = std::vector<unsigned>(param.m_initial_num_cells, 2);
            output = param.ProcessUnsignedVector(input);

            for (unsigned i = 0; i < test_output.size(); i++)
            {
                TS_ASSERT_EQUALS(output.at(i), test_output.at(i));
            }

            // Other vector sizes should raise an error
            input = "2,5,2";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessUnsignedVector(input));
            input = "";
            TS_ASSERT_THROWS_ANYTHING(param.ProcessUnsignedVector(input));
        }

        void TestDeathClockParametersPrint()
        {
            DeathClockParameters param;

            TS_ASSERT_THROWS_NOTHING(std::cout << param.Print(););

            // See if it works as input
            std::string output = param.Print();
            Dictionary dict = create_dictionary(output);;

            DeathClockParameters param2;

            TS_ASSERT_THROWS_NOTHING(param2.Update(dict););

            TS_ASSERT_EQUALS(param2.m_output_directory, "death-clock");
            TS_ASSERT_EQUALS(param2.m_simulation_id, "0");

            TS_ASSERT_EQUALS(param2.m_log_level, 1);

            TS_ASSERT_EQUALS(param2.m_seed, 0);

            TS_ASSERT_EQUALS(param2.m_simulation_time, 100);
            TS_ASSERT_EQUALS(param2.m_dt, 0.01);
            TS_ASSERT_EQUALS(param2.m_sampling_timestep_multiple, 100);

            TS_ASSERT_EQUALS(param2.m_cell_rearrangement_threshold, 0.01);
            TS_ASSERT_EQUALS(param2.m_cell_rearrangement_ratio, 1.5);
            TS_ASSERT_EQUALS(param2.m_t2_threshold, 0.001);

            TS_ASSERT_EQUALS(param2.m_min_cell_count, 0);
            TS_ASSERT_EQUALS(param2.m_max_cell_count, static_cast<unsigned>(-1));

            TS_ASSERT_EQUALS(param2.m_min_cell_count_for_ancestor.size(), 0);
            TS_ASSERT_EQUALS(param2.m_max_cell_count_for_ancestor.size(), 0);

            TS_ASSERT_EQUALS(param2.m_num_cells_across, 4);
            TS_ASSERT_EQUALS(param2.m_num_cells_up, 4);

            TS_ASSERT_EQUALS(param2.m_initial_num_cells, 16);

            std::vector<unsigned> ancestors_test;
            for (unsigned i = 0; i < param2.m_initial_num_cells; i++)
            {
                ancestors_test.push_back(i);
            }
            TS_ASSERT_EQUALS_VECTORS(param2.m_ancestors, ancestors_test);

            TS_ASSERT_EQUALS(param2.m_force_model, "farhadifar");

            TS_ASSERT_EQUALS(param2.m_boundary_tension_parameter, 0.12);

            TS_ASSERT_EQUALS_VECTORS(param2.m_cell_cycle_models,
                    std::vector<unsigned>(param2.m_initial_num_cells,
                        exponential));

            TS_ASSERT_EQUALS_VECTORS(param2.m_reference_target_areas,
                    std::vector<double>(param2.m_initial_num_cells, 1.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_g1_durations,
                    std::vector<double>(param2.m_initial_num_cells, 50.0));
            TS_ASSERT_EQUALS_VECTORS(param2.m_g2_durations,
                    std::vector<double>(param2.m_initial_num_cells, 50.0));
            TS_ASSERT_EQUALS_VECTORS(param2.m_uniform_g1_ranges,
                    std::vector<double>(param2.m_initial_num_cells, 10.0));
            TS_ASSERT_EQUALS_VECTORS(param2.m_max_transit_generations,
                    std::vector<unsigned>(param2.m_initial_num_cells,
                        static_cast<unsigned>(-1)));

            TS_ASSERT_EQUALS_VECTORS(param2.m_apoptosis_times,
                    std::vector<double>(param2.m_initial_num_cells, 0.25));
            TS_ASSERT_EQUALS_VECTORS(param2.m_death_thresholds,
                    std::vector<double>(param2.m_initial_num_cells, 100.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_neighbours_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_neighbours_g2_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));
            TS_ASSERT_EQUALS_VECTORS(param2.m_neighbours_apoptosis_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_normalised_neighbours_g2_constants,
                    std::vector<double>(param2.m_initial_num_cells, 6.0));
            TS_ASSERT_EQUALS_VECTORS(param2.m_normalised_neighbours_apoptosis_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_shared_edge_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_shared_edge_g2_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));
            TS_ASSERT_EQUALS_VECTORS(param2.m_shared_edge_apoptosis_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_normalised_shared_edge_g2_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));
            TS_ASSERT_EQUALS_VECTORS(param2.m_normalised_shared_edge_apoptosis_constants,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_base_rates,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param2.m_initial_taus,
                    std::vector<double>(param2.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS(param2.m_birth_times.size(), param2.m_initial_num_cells);
        }

        void TestDeathClockParametersDefault()
        {
            DeathClockParameters param;

            TS_ASSERT_EQUALS(param.m_output_directory, "death-clock");
            TS_ASSERT_EQUALS(param.m_simulation_id, "0");

            TS_ASSERT_EQUALS(param.m_log_level, 1);

            TS_ASSERT_EQUALS(param.m_seed, 0);

            TS_ASSERT_EQUALS(param.m_simulation_time, 100);
            TS_ASSERT_EQUALS(param.m_dt, 0.01);
            TS_ASSERT_EQUALS(param.m_sampling_timestep_multiple, 100);

            TS_ASSERT_EQUALS(param.m_cell_rearrangement_threshold, 0.01);
            TS_ASSERT_EQUALS(param.m_cell_rearrangement_ratio, 1.5);
            TS_ASSERT_EQUALS(param.m_t2_threshold, 0.001);

            TS_ASSERT_EQUALS(param.m_min_cell_count, 0);
            TS_ASSERT_EQUALS(param.m_max_cell_count, static_cast<unsigned>(-1));

            TS_ASSERT_EQUALS(param.m_min_cell_count_for_ancestor.size(), 0);
            TS_ASSERT_EQUALS(param.m_max_cell_count_for_ancestor.size(), 0);

            TS_ASSERT_EQUALS(param.m_num_cells_across, 4);
            TS_ASSERT_EQUALS(param.m_num_cells_up, 4);

            TS_ASSERT_EQUALS(param.m_initial_num_cells, 16);

            std::vector<unsigned> ancestors_test;
            for (unsigned i = 0; i < param.m_initial_num_cells; i++)
            {
                ancestors_test.push_back(i);
            }
            TS_ASSERT_EQUALS_VECTORS(param.m_ancestors, ancestors_test);

            TS_ASSERT_EQUALS(param.m_force_model, "farhadifar");

            TS_ASSERT_EQUALS(param.m_boundary_tension_parameter, 0.12);

            TS_ASSERT_EQUALS_VECTORS(param.m_cell_cycle_models,
                    std::vector<unsigned>(param.m_initial_num_cells,
                        exponential));

            TS_ASSERT_EQUALS_VECTORS(param.m_reference_target_areas,
                    std::vector<double>(param.m_initial_num_cells, 1.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_g1_durations,
                    std::vector<double>(param.m_initial_num_cells, 50.0));
            TS_ASSERT_EQUALS_VECTORS(param.m_g2_durations,
                    std::vector<double>(param.m_initial_num_cells, 50.0));
            TS_ASSERT_EQUALS_VECTORS(param.m_uniform_g1_ranges,
                    std::vector<double>(param.m_initial_num_cells, 10.0));
            TS_ASSERT_EQUALS_VECTORS(param.m_max_transit_generations,
                    std::vector<unsigned>(param.m_initial_num_cells,
                        static_cast<unsigned>(-1)));

            TS_ASSERT_EQUALS_VECTORS(param.m_apoptosis_times,
                    std::vector<double>(param.m_initial_num_cells, 0.25));
            TS_ASSERT_EQUALS_VECTORS(param.m_death_thresholds,
                    std::vector<double>(param.m_initial_num_cells, 100.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_neighbours_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_neighbours_g2_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));
            TS_ASSERT_EQUALS_VECTORS(param.m_neighbours_apoptosis_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_normalised_neighbours_g2_constants,
                    std::vector<double>(param.m_initial_num_cells, 6.0));
            TS_ASSERT_EQUALS_VECTORS(param.m_normalised_neighbours_apoptosis_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_shared_edge_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_shared_edge_g2_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));
            TS_ASSERT_EQUALS_VECTORS(param.m_shared_edge_apoptosis_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_normalised_shared_edge_g2_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));
            TS_ASSERT_EQUALS_VECTORS(param.m_normalised_shared_edge_apoptosis_constants,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_base_rates,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS_VECTORS(param.m_initial_taus,
                    std::vector<double>(param.m_initial_num_cells, 0.0));

            TS_ASSERT_EQUALS(param.m_birth_times.size(), param.m_initial_num_cells);
        }

        // Test Update of
        // output-directory
        void TestDeathClockParametersUpdate1()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("output-directory = other-dir");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_output_directory, "other-dir");
        }

        // seed
        void TestDeathClockParametersUpdate2()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("seed = 3");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_seed, 3);
        }

        // simulation-time
        void TestDeathClockParametersUpdate3()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("simulation-time = 50");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_simulation_time, 50.0);
        }

        // num-cells-across
        void TestDeathClockParametersUpdate4()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("num-cells-across = 2");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_num_cells_across, 2);
            TS_ASSERT_EQUALS(param.m_ancestors.size(), 8);
            TS_ASSERT_EQUALS(param.m_cell_cycle_models.size(), 8);
            TS_ASSERT_EQUALS(param.m_reference_target_areas.size(), 8);
            TS_ASSERT_EQUALS(param.m_birth_times.size(), 8);
        }

        // num-cells-up
        void TestDeathClockParametersUpdate5()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("num-cells-up = 2");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_num_cells_up, 2);
            TS_ASSERT_EQUALS(param.m_ancestors.size(), 8);
            TS_ASSERT_EQUALS(param.m_cell_cycle_models.size(), 8);
            TS_ASSERT_EQUALS(param.m_reference_target_areas.size(), 8);
            TS_ASSERT_EQUALS(param.m_birth_times.size(), 8);
        }

        // num-cells-across and num-cells-up
        void TestDeathClockParametersUpdate6()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("num-cells-across = 2\nnum-cells-up = 2");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_num_cells_across, 2);
            TS_ASSERT_EQUALS(param.m_num_cells_up, 2);
            TS_ASSERT_EQUALS(param.m_ancestors.size(), 4);
            TS_ASSERT_EQUALS(param.m_cell_cycle_models.size(), 4);
            TS_ASSERT_EQUALS(param.m_reference_target_areas.size(), 4);
            TS_ASSERT_EQUALS(param.m_birth_times.size(), 4);
        }

        // ancestors
        void TestDeathClockParametersUpdate7()
        {
            DeathClockParameters param;

            // All entries
            Dictionary dict = create_dictionary("ancestors = "
            "0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150");
            param.Update(dict);

            std::vector<unsigned> ancestors_test;
            for (unsigned i = 0; i < param.m_initial_num_cells; i++)
            {
                ancestors_test.push_back(i * 10);
            }
            TS_ASSERT_EQUALS_VECTORS(param.m_ancestors, ancestors_test);

            // Single entry
            dict = create_dictionary("ancestors = 0");
            param.Update(dict);

            TS_ASSERT_EQUALS_VECTORS(param.m_ancestors, std::vector<unsigned>(param.m_initial_num_cells, 0));

            // All entries
            dict = create_dictionary("ancestors = 0, 10, 20, 30\n"
                "num-cells-across = 2\n"
                "num-cells-up = 2");
            param.Update(dict);

            ancestors_test.clear();
            for (unsigned i = 0; i < param.m_initial_num_cells; i++)
            {
                ancestors_test.push_back(i * 10);
            }
            TS_ASSERT_EQUALS_VECTORS(param.m_ancestors, ancestors_test);

            // Single entry
            dict = create_dictionary("ancestors = 0");
            param.Update(dict);

            TS_ASSERT_EQUALS_VECTORS(param.m_ancestors, std::vector<unsigned>(param.m_initial_num_cells, 0));
        }

        // cell-cycle-models
        void TestDeathClockParametersUpdate8()
        {
            DeathClockParameters param;

            // All entries
            Dictionary dict = create_dictionary("cell-cycle-models = "
            "exponential, "
            "exponential, "
            "exponential, "
            "exponential, "
            "exponential, "
            "exponential, "
            "exponential, "
            "exponential, "
            "uniform, "
            "uniform, "
            "uniform, "
            "uniform, "
            "uniform, "
            "uniform, "
            "uniform, "
            "uniform");
            param.Update(dict);

            std::vector<unsigned> cell_cycle_models_test;
            for (unsigned i = 0; i < param.m_initial_num_cells / 2; i++)
            {
                cell_cycle_models_test.push_back(exponential);
            }

            for (unsigned i = 0; i < param.m_initial_num_cells / 2; i++)
            {
                cell_cycle_models_test.push_back(uniform);
            }
            TS_ASSERT_EQUALS_VECTORS(param.m_cell_cycle_models, cell_cycle_models_test);

            // Single entry
            dict = create_dictionary("cell-cycle-models = exponential");
            param.Update(dict);

            TS_ASSERT_EQUALS_VECTORS(param.m_cell_cycle_models, std::vector<unsigned>(param.m_initial_num_cells, exponential));

            // All entries
            dict = create_dictionary("cell-cycle-models = "
            "exponential, "
            "exponential, "
            "uniform, "
            "uniform\n"
            "num-cells-across = 2\n"
            "num-cells-up = 2");
            param.Update(dict);

            cell_cycle_models_test.clear();
            for (unsigned i = 0; i < param.m_initial_num_cells / 2; i++)
            {
                cell_cycle_models_test.push_back(exponential);
            }

            for (unsigned i = 0; i < param.m_initial_num_cells / 2; i++)
            {
                cell_cycle_models_test.push_back(uniform);
            }
            TS_ASSERT_EQUALS_VECTORS(param.m_cell_cycle_models, cell_cycle_models_test);

            // Single entry
            dict = create_dictionary("cell-cycle-models = exponential");
            param.Update(dict);

            TS_ASSERT_EQUALS_VECTORS(param.m_cell_cycle_models, std::vector<unsigned>(param.m_initial_num_cells, exponential));
        }

        // reference-target-areas
        void TestDeathClockParametersUpdate9()
        {
            DeathClockParameters param;

            // All entries
            Dictionary dict = create_dictionary("reference-target-areas = "
            "0.3, 10.3, 20.3, 30.3, 40.3, 50.3, 60.3, 70.3, 80.3, 90.3, 100.3, 110.3, 120.3, 130.3, 140.3, 150.3");
            param.Update(dict);

            std::vector<double> reference_target_areas_test;
            for (unsigned i = 0; i < param.m_initial_num_cells; i++)
            {
                reference_target_areas_test.push_back(static_cast<double>(i * 10) + 0.3);
            }
            TS_ASSERT_EQUALS_VECTORS(param.m_reference_target_areas, reference_target_areas_test);

            // Single entry
            dict = create_dictionary("reference-target-areas = 0.3");
            param.Update(dict);

            TS_ASSERT_EQUALS_VECTORS(param.m_reference_target_areas, std::vector<double>(param.m_initial_num_cells, 0.3));

            // All entries
            dict = create_dictionary("reference-target-areas = "
            "0.3, 10.3, 20.3, 30.3\n"
                "num-cells-across = 2\n"
                "num-cells-up = 2");
            param.Update(dict);

            reference_target_areas_test.clear();
            for (unsigned i = 0; i < param.m_initial_num_cells; i++)
            {
                reference_target_areas_test.push_back(static_cast<double>(i * 10) + 0.3);
            }
            TS_ASSERT_EQUALS_VECTORS(param.m_reference_target_areas, reference_target_areas_test);

            // Single entry
            dict = create_dictionary("reference-target-areas = 0.3");
            param.Update(dict);

            TS_ASSERT_EQUALS_VECTORS(param.m_reference_target_areas, std::vector<double>(param.m_initial_num_cells, 0.3));
        }

        // g1-durations and g2-durations
        void TestDeathClockParametersUpdate10()
        {
            DeathClockParameters param;

            // All entries
            Dictionary dict = create_dictionary("g1-durations = "
            "0.3, 10.3, 20.3, 30.3, 40.3, 50.3, 60.3, 70.3, 80.3, 90.3, 100.3, 110.3, 120.3, 130.3, 140.3, 150.3\n"
            "g2-durations = "
            "0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0");

            param.Update(dict);

            std::vector<double> birth_time_mins;
            for (unsigned i = 0; i < param.m_initial_num_cells; i++)
            {
                birth_time_mins.push_back(-(static_cast<double>(i * 11) + 0.3));
            }
            TS_ASSERT_LESS_THAN_EQUALS_VECTORS(birth_time_mins, param.m_birth_times);

            // Single entry
            dict = create_dictionary("g1-durations = 0.3\ng2-durations = 10.0");
            param.Update(dict);

            TS_ASSERT_LESS_THAN_EQUALS_VECTORS(std::vector<double>(param.m_initial_num_cells,
                        -10.3), param.m_birth_times);

            // All entries
            dict = create_dictionary("g1-durations = "
            "0.3, 10.3, 20.3, 30.3\n"
            "g2-durations = "
            "0.0, 1.0, 2.0, 3.0\n"
            "num-cells-across = 2\n"
            "num-cells-up = 2");

            param.Update(dict);

            birth_time_mins.clear();
            for (unsigned i = 0; i < param.m_initial_num_cells; i++)
            {
                birth_time_mins.push_back(-(static_cast<double>(i * 11) + 0.3));
            }
            TS_ASSERT_LESS_THAN_EQUALS_VECTORS(birth_time_mins, param.m_birth_times);

            // Single entry
            dict = create_dictionary("g1-durations = 0.3\ng2-durations = 10.0");
            param.Update(dict);

            TS_ASSERT_LESS_THAN_EQUALS_VECTORS(std::vector<double>(param.m_initial_num_cells,
                        -10.3), param.m_birth_times);
        }

        // force-model
        void TestDeathClockParametersUpdate11()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("force-model = nagai-honda");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_force_model, "nagai-honda");
        }

        // min-cell-count-for-ancestor and max-cell-count-for-ancestor
        void TestDeathClockParametersUpdate12()
        {
            DeathClockParameters param;

            Dictionary dict = create_dictionary("min-cell-count-for-ancestor = "
                    "{0:20, 323: 230}\n"
                    "max-cell-count-for-ancestor = "
                    "{23: 3221, 4343: 342342}");
            param.Update(dict);

            TS_ASSERT_EQUALS(param.m_min_cell_count_for_ancestor[0], 20);
            TS_ASSERT_EQUALS(param.m_min_cell_count_for_ancestor[323], 230);
            TS_ASSERT_EQUALS(param.m_max_cell_count_for_ancestor[23], 3221);
            TS_ASSERT_EQUALS(param.m_max_cell_count_for_ancestor[4343], 342342);
        }

        // wrong cell cycle model
        void TestDeathClockParametersError1()
        {
            DeathClockParameters param;
            Dictionary dict = create_dictionary("cell-cycle-models = "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););

            dict = create_dictionary("cell-cycle-models = uniform, exponential");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););
        }

        // incorrect unsigned
        void TestDeathClockParametersError2()
        {
            DeathClockParameters param;
            Dictionary dict = create_dictionary("seed = rofl");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););
        }

        // incorrect double
        void TestDeathClockParametersError3()
        {
            DeathClockParameters param;
            Dictionary dict = create_dictionary("simulation-time = lol");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););
        }

        // incorrect vector of unsigneds
        void TestDeathClockParametersError4()
        {
            DeathClockParameters param;
            Dictionary dict = create_dictionary("ancestors = "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););

            dict = create_dictionary("ancestors = 0, 1");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););
        }

        // incorrect vector of doubles
        void TestDeathClockParametersError5()
        {
            DeathClockParameters param;
            Dictionary dict = create_dictionary("reference-target-areas = "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao, "
            "lmao");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););

            dict = create_dictionary("reference-target-areas = 2.0, 1.5");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););
        }

        // Incorrect unsigned to unsigned map
        void TestDeathClockParametersError6()
        {
            DeathClockParameters param;
            Dictionary dict = create_dictionary("min-cell-count-for-ancestor = "
                    "132 : 5343, 4343 : 535432");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););

            dict = create_dictionary("min-cell-count-for-ancestor = "
                    "{lmao : haha, hihi : hoho}");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););

            dict = create_dictionary("min-cell-count-for-ancestor = "
                    "{0 1, 1: 2}");
            TS_ASSERT_THROWS_ANYTHING(param.Update(dict););
        }

        void TestDeathClockParametersHelp()
        {
            TS_ASSERT_THROWS_NOTHING(std::cout << DeathClockParameters::Help(););
        }

        void TestDeathClockParametersIsValid()
        {
            DeathClockParameters param;

            // Default state should be valid
            TS_ASSERT_EQUALS(param.IsValid(), true);

            // Changing m_initial_num_cells should invalidate it
            param.m_initial_num_cells = 2;
            TS_ASSERT_EQUALS(param.IsValid(), false);

            // Update with an empty dictionary should fix this
            Dictionary dict = create_dictionary("");
            param.Update(dict);
            TS_ASSERT_EQUALS(param.IsValid(), true);

            // Changing m_initial_num_cells, then updating with non-empty
            // dictionary should also leave it valid
            param.m_initial_num_cells = 2;
            dict = create_dictionary("g1-duration = 30");
            param.Update(dict);
            TS_ASSERT_EQUALS(param.IsValid(), true);

            // You can't trick the struct into thinking no change has occured
            param.m_initial_num_cells = 2;
            dict = create_dictionary("num-cells-across=2\nnum-cells-up=1");
            param.Update(dict);
            TS_ASSERT_EQUALS(param.IsValid(), true);

            // Changing the size of a vector should invalid the struct
            param.m_ancestors.pop_back();
            TS_ASSERT_EQUALS(param.IsValid(), false);

            // Changing a birth time so that it is positive or further in the
            // past than the average g1, g2 durations combined will also
            // invalid the struct
            param.Update(dict);
            TS_ASSERT_EQUALS(param.IsValid(), true);
            param.m_birth_times[0] = 1.0;
            TS_ASSERT_EQUALS(param.IsValid(), false);
            param.Update(dict);
            TS_ASSERT_EQUALS(param.IsValid(), true);
            param.m_birth_times[0] = - 2.0 * (param.m_g1_durations[0] + param.m_g2_durations[0]);
            TS_ASSERT_EQUALS(param.IsValid(), false);
        }
};

#endif // TESTDEATHCLOCKPARAMETERS_HPP
