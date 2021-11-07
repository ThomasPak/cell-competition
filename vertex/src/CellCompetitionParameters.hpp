#ifndef CELLCOMPETITIONPARAMETERS_HPP_
#define CELLCOMPETITIONPARAMETERS_HPP_

#include <string>

#include "Dictionary.hpp"

/** Struct that contains parameters for CellCompetitionSimulation().
 *
 * @todo Document every member.
 * @todo Rename every member to conform with PascalCase.
 */
struct CellCompetitionParameters
{
    /** Output directory.  Defaults to "cell-competition". */
    std::string m_output_directory = "cell-competition";
    std::string m_simulation_id = "0";

    unsigned m_seed = 0;

    double m_simulation_time = 100;
    double m_dt = 0.01;
    unsigned m_sampling_timestep_multiple = 100;

    double m_cell_rearrangement_threshold = 0.01;
    double m_cell_rearrangement_ratio = 1.5;
    double m_t2_threshold = 0.001;

    unsigned m_num_cells_across = 4;
    unsigned m_num_cells_up = 4;
    bool m_random_labelling = true;
    double m_initial_b_ratio = 0.5;

    double m_random_movement_parameter = 0.1;
    double m_cell_ab_line_tension_parameter = 0.12;

    // Cell A
    double m_cell_a_target_area = 1.0;
    double m_cell_a_elasticity_parameter = 1.0;
    double m_cell_a_contractility_parameter = 0.04;
    double m_cell_a_line_tension_parameter = 0.12;
    double m_cell_a_boundary_tension_parameter = 0.12;
    double m_cell_a_g1_duration = 30.0;
    double m_cell_a_g2_duration = 70.0;
    unsigned m_cell_a_max_generations = static_cast<unsigned>(-1);

    // Cell B
    double m_cell_b_target_area = 1.0;
    double m_cell_b_elasticity_parameter = 1.0;
    double m_cell_b_contractility_parameter = 0.04;
    double m_cell_b_line_tension_parameter = 0.12;
    double m_cell_b_boundary_tension_parameter = 0.12;
    double m_cell_b_g1_duration = 30.0;
    double m_cell_b_g2_duration = 70.0;
    unsigned m_cell_b_max_generations = static_cast<unsigned>(-1);

    /** Update values given by a Dictionary.
     *
     * @param dict  Dictionary containing updated values.
     */
    void Update(const Dictionary& dict);

    /** Print parameter values.
     *
     * @return string containing printed parameter values.
     */
    std::string Print();

    /** Generate Help message
     *
     * @return help string.
     */
    static std::string Help();
};

#endif // CELLCOMPETITIONPARAMETERS_HPP_
