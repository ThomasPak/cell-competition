#include <sstream>

#include "Dictionary.hpp"
#include "CellCompetitionParameters.hpp"

void CellCompetitionParameters::Update(const Dictionary& dict)
{
    if (dict.count("output-directory") == 1)
        m_output_directory = dict.at("output-directory");

    if (dict.count("simulation-id") == 1)
        m_simulation_id = dict.at("simulation-id");

    if (dict.count("seed") == 1)
        m_seed = std::stoul(dict.at("seed"));

    if (dict.count("simulation-time") == 1)
        m_simulation_time = std::stod(dict.at("simulation-time"));

    if (dict.count("dt") == 1)
        m_dt = std::stod(dict.at("dt"));

    if (dict.count("sampling-timestep-multiple") == 1)
        m_sampling_timestep_multiple =
            std::stoul(dict.at("sampling-timestep-multiple"));

    if (dict.count("cell-rearrangement-threshold") == 1)
        m_cell_rearrangement_threshold =
            std::stod(dict.at("cell-rearrangement-threshold"));

    if (dict.count("cell-rearrangement-ratio") == 1)
        m_cell_rearrangement_ratio =
            std::stod(dict.at("cell-rearrangement-ratio"));

    if (dict.count("t2-threshold") == 1)
        m_t2_threshold =
            std::stod(dict.at("t2-threshold"));

    if (dict.count("num-cells-across") == 1)
        m_num_cells_across = std::stoul(dict.at("num-cells-across"));

    if (dict.count("num-cells-up") == 1)
        m_num_cells_up = std::stoul(dict.at("num-cells-up"));

    if (dict.count("random-movement-parameter") == 1)
        m_random_movement_parameter = std::stod(dict.at("random-movement-parameter"));

    if (dict.count("cell-ab-line-tension-parameter") == 1)
        m_cell_ab_line_tension_parameter =
            std::stod(dict.at("cell-ab-line-tension-parameter"));

    if (dict.count("random-labelling") == 1)
    {
        if (dict.at("random-labelling") == "true")
            m_random_labelling = true;
        else if (dict.at("random-labelling") == "false")
            m_random_labelling = false;
        else
        {
            std::string err_msg = "Error: value of random-labelling must be `true` or `false`";
            throw std::runtime_error(err_msg);
        }
    }

    if (dict.count("initial-b-ratio") == 1)
        m_initial_b_ratio = std::stod(dict.at("initial-b-ratio"));

    // Cell type A
    if (dict.count("cell-a-target-area") == 1)
        m_cell_a_target_area = std::stod(dict.at("cell-a-target-area"));

    if (dict.count("cell-a-elasticity-parameter") == 1)
        m_cell_a_elasticity_parameter =
            std::stod(dict.at("cell-a-elasticity-parameter"));

    if (dict.count("cell-a-contractility-parameter") == 1)
        m_cell_a_contractility_parameter =
            std::stod(dict.at("cell-a-contractility-parameter"));

    if (dict.count("cell-a-line-tension-parameter") == 1)
        m_cell_a_line_tension_parameter =
            std::stod(dict.at("cell-a-line-tension-parameter"));

    if (dict.count("cell-a-boundary-tension-parameter") == 1)
        m_cell_a_boundary_tension_parameter =
            std::stod(dict.at("cell-a-boundary-tension-parameter"));

    if (dict.count("cell-a-g1-duration") == 1)
        m_cell_a_g1_duration = std::stod(dict.at("cell-a-g1-duration"));

    if (dict.count("cell-a-g2-duration") == 1)
        m_cell_a_g2_duration = std::stod(dict.at("cell-a-g2-duration"));

    if (dict.count("cell-a-max-generations") == 1)
        m_cell_a_max_generations =
            std::stoul(dict.at("cell-a-max-generations"));

    // Cell type B
    if (dict.count("cell-b-target-area") == 1)
        m_cell_b_target_area = std::stod(dict.at("cell-b-target-area"));

    if (dict.count("cell-b-elasticity-parameter") == 1)
        m_cell_b_elasticity_parameter =
            std::stod(dict.at("cell-b-elasticity-parameter"));

    if (dict.count("cell-b-contractility-parameter") == 1)
        m_cell_b_contractility_parameter =
            std::stod(dict.at("cell-b-contractility-parameter"));

    if (dict.count("cell-b-line-tension-parameter") == 1)
        m_cell_b_line_tension_parameter =
            std::stod(dict.at("cell-b-line-tension-parameter"));

    if (dict.count("cell-b-boundary-tension-parameter") == 1)
        m_cell_b_boundary_tension_parameter =
            std::stod(dict.at("cell-b-boundary-tension-parameter"));

    if (dict.count("cell-b-g1-duration") == 1)
        m_cell_b_g1_duration = std::stod(dict.at("cell-b-g1-duration"));

    if (dict.count("cell-b-g2-duration") == 1)
        m_cell_b_g2_duration = std::stod(dict.at("cell-b-g2-duration"));

    if (dict.count("cell-b-max-generations") == 1)
        m_cell_b_max_generations =
            std::stoul(dict.at("cell-b-max-generations"));
}

std::string CellCompetitionParameters::Print()
{
    std::stringstream print_sstrm;
    print_sstrm << "Cell competition parameters:\n";

    print_sstrm << std::endl;
    print_sstrm << "output-directory: " << m_output_directory << std::endl;
    print_sstrm << "simulation-id: " << m_simulation_id << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "seed: " << m_seed << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "simulation-time: " << m_simulation_time << std::endl;
    print_sstrm << "dt: " << m_dt << std::endl;
    print_sstrm << "sampling-timestep-multiple: " <<
        m_sampling_timestep_multiple << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "cell-rearrangement-threshold: " <<
        m_cell_rearrangement_threshold << std::endl;
    print_sstrm << "cell-rearrangement-ratio: " << m_cell_rearrangement_ratio
        << std::endl;
    print_sstrm << "t2-threshold: " << m_t2_threshold << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "num-cells-across: " << m_num_cells_across << std::endl;
    print_sstrm << "num-cells-up: " << m_num_cells_up << std::endl;
    print_sstrm << "random-labelling: " << (m_random_labelling == true ? "true" :
            "false") << std::endl;
    print_sstrm << "initial-b-ratio: " << m_initial_b_ratio << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "random-movement-parameter: " << m_random_movement_parameter
        << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "cell-ab-line-tension-parameter: " <<
        m_cell_ab_line_tension_parameter << std::endl;

    // Cell A
    print_sstrm << std::endl;
    print_sstrm << "cell-a-target-area: " << m_cell_a_target_area << std::endl;
    print_sstrm << "cell-a-elasticity-parameter: " <<
        m_cell_a_elasticity_parameter << std::endl;
    print_sstrm << "cell-a-contractility-parameter: " <<
        m_cell_a_contractility_parameter << std::endl;
    print_sstrm << "cell-a-line-tension-parameter: " <<
        m_cell_a_line_tension_parameter << std::endl;
    print_sstrm << "cell-a-boundary-tension-parameter: " <<
        m_cell_a_boundary_tension_parameter << std::endl;
    print_sstrm << "cell-a-g1-duration: " << m_cell_a_g1_duration << std::endl;
    print_sstrm << "cell-a-g2-duration: " << m_cell_a_g2_duration << std::endl;
    print_sstrm << "cell-a-max-generations: " << m_cell_a_max_generations <<
        std::endl;

    // Cell B
    print_sstrm << std::endl;
    print_sstrm << "cell-b-target-area: " << m_cell_b_target_area << std::endl;
    print_sstrm << "cell-b-elasticity-parameter: " <<
        m_cell_b_elasticity_parameter << std::endl;
    print_sstrm << "cell-b-contractility-parameter: " <<
        m_cell_b_contractility_parameter << std::endl;
    print_sstrm << "cell-b-line-tension-parameter: " <<
        m_cell_b_line_tension_parameter << std::endl;
    print_sstrm << "cell-b-boundary-tension-parameter: " <<
        m_cell_b_boundary_tension_parameter << std::endl;
    print_sstrm << "cell-b-g1-duration: " << m_cell_b_g1_duration << std::endl;
    print_sstrm << "cell-b-g2-duration: " << m_cell_b_g2_duration << std::endl;
    print_sstrm << "cell-b-max-generations: " << m_cell_b_max_generations <<
        std::endl;

    return print_sstrm.str();
}

std::string CellCompetitionParameters::Help()
{
    std::string help_msg("Valid keys:\n");

    help_msg += "\n";
    help_msg += "  output-directory                     (default: cell-competition)\n";
    help_msg += "  simulation-id                        (default: 0)\n";

    help_msg += "\n";
    help_msg += "  seed                                 (default: 0 (-> random))\n";

    help_msg += "\n";
    help_msg += "  simulation-time                      (default: 100)\n";
    help_msg += "  dt                                   (default: 0.01)\n";
    help_msg += "  sampling-timestep-multiple           (default: 100)\n";

    help_msg += "\n";
    help_msg += "  cell-rearrangement-threshold         (default: 0.01)\n";
    help_msg += "  cell-rearrangement-ratio             (default: 1.5)\n";
    help_msg += "  t2-threshold                         (default: 0.001)\n";

    help_msg += "\n";
    help_msg += "  num-cells-across                     (default: 4)\n";
    help_msg += "  num-cells-up                         (default: 4)\n";
    help_msg += "  random-labelling (true|false)        (default: true)\n";
    help_msg += "  initial-b-ratio                      (default: 0.5)\n";

    help_msg += "\n";
    help_msg += "  random-movement-parameter            (default: 0.1)\n";

    help_msg += "\n";
    help_msg += "  cell-ab-line-tension-parameter       (default: 0.12)\n";

    // Cell A
    help_msg += "\n";
    help_msg += "  cell-a-target-area                   (default: 1.0)\n";
    help_msg += "  cell-a-elasticity-parameter          (default: 1.0)\n";
    help_msg += "  cell-a-contractility-parameter       (default: 0.04)\n";
    help_msg += "  cell-a-line-tension-parameter        (default: 0.12)\n";
    help_msg += "  cell-a-boundary-tension-parameter    (default: 0.12)\n";
    help_msg += "  cell-a-g1-duration                   (default: 30.0)\n";
    help_msg += "  cell-a-g2-duration                   (default: 70.0)\n";
    help_msg += "  cell-a-max-generations               (default: 3)\n";

    // Cell B
    help_msg += "\n";
    help_msg += "  cell-b-target-area                   (default: 1.0)\n";
    help_msg += "  cell-b-elasticity-parameter          (default: 1.0)\n";
    help_msg += "  cell-b-contractility-parameter       (default: 0.04)\n";
    help_msg += "  cell-b-line-tension-parameter        (default: 0.12)\n";
    help_msg += "  cell-b-boundary-tension-parameter    (default: 0.12)\n";
    help_msg += "  cell-b-g1-duration                   (default: 30.0)\n";
    help_msg += "  cell-b-g2-duration                   (default: 70.0)\n";
    help_msg += "  cell-b-max-generations               (default: 3)\n";

    return help_msg;
}
