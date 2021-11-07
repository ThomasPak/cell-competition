#include <sstream>
#include <string>
#include <map>

#include "Exception.hpp"
#include "RandomNumberGenerator.hpp"

#include "Dictionary.hpp"
#include "DeathClockParameters.hpp"

CellCycleModel StringToCellCycleModel(const std::string& input)
{
    if (input.compare("exponential") == 0)
	{
        return exponential;
	}
    else if (input.compare("uniform") == 0)
	{
        return uniform;
	}

	std::string error_msg;
	error_msg += "\"";
	error_msg += input;
	error_msg += "\" could not be interpreted as a cell cycle model";
	EXCEPTION(error_msg);
}

std::string CellCycleModelToString(CellCycleModel ccm)
{
	switch (ccm)
	{
		case exponential:
			return "exponential";
		case uniform:
			return "uniform";
	}

	NEVER_REACHED;
}

DeathClockParameters::DeathClockParameters()
{
    // Initialise ancestors
    for (unsigned ancestor = 0; ancestor < m_initial_num_cells; ancestor++)
    {
        m_ancestors.push_back(ancestor);
    }

    // Initialise birth times
    for (unsigned i = 0; i < m_initial_num_cells; i++)
    {
        double average_cc_duration = m_g1_durations[i] + m_g2_durations[i];
        double birth_time = - average_cc_duration *
            RandomNumberGenerator::Instance()->ranf();
        m_birth_times.push_back(birth_time);
    }
}

void DeathClockParameters::Update(const Dictionary& dict)
{
    if (dict.count("output-directory") == 1)
        m_output_directory = dict.at("output-directory");

    if (dict.count("simulation-id") == 1)
        m_simulation_id = dict.at("simulation-id");

    if (dict.count("log-level") == 1)
        m_log_level = std::stoul(dict.at("log-level"));

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

    if (dict.count("min-cell-count") == 1)
        m_min_cell_count = std::stoul(dict.at("min-cell-count"));

    if (dict.count("max-cell-count") == 1)
        m_max_cell_count = std::stoul(dict.at("max-cell-count"));

    if (dict.count("min-cell-count-for-ancestor") == 1)
        m_min_cell_count_for_ancestor =
            ReadUnsignedToUnsignedMap(dict.at("min-cell-count-for-ancestor"));

    if (dict.count("max-cell-count-for-ancestor") == 1)
        m_max_cell_count_for_ancestor =
            ReadUnsignedToUnsignedMap(dict.at("max-cell-count-for-ancestor"));

    if (dict.count("num-cells-across") == 1)
        m_num_cells_across = std::stoul(dict.at("num-cells-across"));

    if (dict.count("num-cells-up") == 1)
        m_num_cells_up = std::stoul(dict.at("num-cells-up"));

    if (dict.count("force-model") == 1)
        m_force_model = dict.at("force-model");

    if (dict.count("boundary-tension-parameter") == 1)
        m_boundary_tension_parameter = std::stod(dict.at("boundary-tension-parameter"));

    // If the struct is no longer valid, reinitialise vectors to restore
    // validity.
    if (!IsValid())
    {
        m_initial_num_cells = m_num_cells_across * m_num_cells_up;

        // Initialise ancestors
        m_ancestors.clear();
        for (unsigned ancestor = 0; ancestor < m_initial_num_cells; ancestor++)
        {
            m_ancestors.push_back(ancestor);
        }

        m_cell_cycle_models = std::vector<unsigned>(m_initial_num_cells, m_cell_cycle_models[0]);

        m_reference_target_areas = std::vector<double>(m_initial_num_cells, m_reference_target_areas[0]);

        m_g1_durations = std::vector<double>(m_initial_num_cells, m_g1_durations[0]);
        m_g2_durations = std::vector<double>(m_initial_num_cells, m_g2_durations[0]);
        m_uniform_g1_ranges = std::vector<double>(m_initial_num_cells, m_uniform_g1_ranges[0]);
        m_max_transit_generations = std::vector<unsigned>(m_initial_num_cells, m_max_transit_generations[0]);

        m_apoptosis_times = std::vector<double>(m_initial_num_cells, m_apoptosis_times[0]);
        m_death_thresholds = std::vector<double>(m_initial_num_cells, m_death_thresholds[0]);

        m_neighbours_constants = std::vector<double>(m_initial_num_cells, m_neighbours_constants[0]);

        m_neighbours_g2_constants = std::vector<double>(m_initial_num_cells, m_neighbours_g2_constants[0]);
        m_neighbours_apoptosis_constants = std::vector<double>(m_initial_num_cells, m_neighbours_apoptosis_constants[0]);

        m_normalised_neighbours_g2_constants = std::vector<double>(m_initial_num_cells, m_normalised_neighbours_g2_constants[0]);
        m_normalised_neighbours_apoptosis_constants = std::vector<double>(m_initial_num_cells, m_normalised_neighbours_apoptosis_constants[0]);

        m_shared_edge_constants = std::vector<double>(m_initial_num_cells, m_shared_edge_constants[0]);

        m_shared_edge_g2_constants = std::vector<double>(m_initial_num_cells, m_shared_edge_g2_constants[0]);
        m_shared_edge_apoptosis_constants = std::vector<double>(m_initial_num_cells, m_shared_edge_apoptosis_constants[0]);

        m_normalised_shared_edge_g2_constants = std::vector<double>(m_initial_num_cells, m_normalised_shared_edge_g2_constants[0]);
        m_normalised_shared_edge_apoptosis_constants = std::vector<double>(m_initial_num_cells, m_normalised_shared_edge_apoptosis_constants[0]);

        m_base_rates = std::vector<double>(m_initial_num_cells, m_base_rates[0]);

        // Define initial state
        m_initial_taus = std::vector<double>(m_initial_num_cells, m_initial_taus[0]);

        // Initialise birth times
        m_birth_times.clear();
        for (unsigned i = 0; i < m_initial_num_cells; i++)
        {
            double average_cc_duration = m_g1_durations[i] + m_g2_durations[i];
            double birth_time = - average_cc_duration *
                RandomNumberGenerator::Instance()->ranf();
            m_birth_times.push_back(birth_time);
        }
    }

    if (dict.count("ancestors") == 1)
        m_ancestors = ProcessUnsignedVector(dict.at("ancestors"));

    if (dict.count("cell-cycle-models") == 1)
        m_cell_cycle_models =
            ProcessUnsignedVector(FormatCellCycleModelToUnsigned(dict.at("cell-cycle-models")));

    if (dict.count("reference-target-areas") == 1)
        m_reference_target_areas = ProcessDoubleVector(dict.at("reference-target-areas"));

    if (dict.count("g1-durations") == 1)
        m_g1_durations = ProcessDoubleVector(dict.at("g1-durations"));

    if (dict.count("g2-durations") == 1)
        m_g2_durations = ProcessDoubleVector(dict.at("g2-durations"));

    // If g1-durations or g2-durations was altered, reinitialise birth times
    if ((dict.count("g1-durations") == 1) || (dict.count("g2-durations") == 1))
    {
        m_birth_times.clear();
        for (unsigned i = 0; i < m_initial_num_cells; i++)
        {
            double average_cc_duration = m_g1_durations[i] + m_g2_durations[i];
            double birth_time = - average_cc_duration *
                RandomNumberGenerator::Instance()->ranf();
            m_birth_times.push_back(birth_time);
        }
    }

    if (dict.count("uniform-g1-ranges") == 1)
        m_uniform_g1_ranges = ProcessDoubleVector(dict.at("uniform-g1-ranges"));

    if (dict.count("max-transit-generations") == 1)
        m_max_transit_generations = ProcessUnsignedVector(dict.at("max-transit-generations"));

    if (dict.count("apoptosis-times") == 1)
        m_apoptosis_times = ProcessDoubleVector(dict.at("apoptosis-times"));

    if (dict.count("death-thresholds") == 1)
        m_death_thresholds = ProcessDoubleVector(dict.at("death-thresholds"));

    if (dict.count("neighbours-constants") == 1)
        m_neighbours_constants = ProcessDoubleVector(dict.at("neighbours-constants"));

    if (dict.count("neighbours-g2-constants") == 1)
        m_neighbours_g2_constants = ProcessDoubleVector(dict.at("neighbours-g2-constants"));

    if (dict.count("neighbours-apoptosis-constants") == 1)
        m_neighbours_apoptosis_constants = ProcessDoubleVector(dict.at("neighbours-apoptosis-constants"));

    if (dict.count("normalised-neighbours-g2-constants") == 1)
        m_normalised_neighbours_g2_constants = ProcessDoubleVector(dict.at("normalised-neighbours-g2-constants"));

    if (dict.count("normalised-neighbours-apoptosis-constants") == 1)
        m_normalised_neighbours_apoptosis_constants = ProcessDoubleVector(dict.at("normalised-neighbours-apoptosis-constants"));

    if (dict.count("shared-edge-constants") == 1)
        m_shared_edge_constants = ProcessDoubleVector(dict.at("shared-edge-constants"));

    if (dict.count("shared-edge-g2-constants") == 1)
        m_shared_edge_g2_constants = ProcessDoubleVector(dict.at("shared-edge-g2-constants"));

    if (dict.count("shared-edge-apoptosis-constants") == 1)
        m_shared_edge_apoptosis_constants = ProcessDoubleVector(dict.at("shared-edge-apoptosis-constants"));

    if (dict.count("normalised-shared-edge-g2-constants") == 1)
        m_normalised_shared_edge_g2_constants = ProcessDoubleVector(dict.at("normalised-shared-edge-g2-constants"));

    if (dict.count("normalised-shared-edge-apoptosis-constants") == 1)
        m_normalised_shared_edge_apoptosis_constants = ProcessDoubleVector(dict.at("normalised-shared-edge-apoptosis-constants"));

    if (dict.count("base-rates") == 1)
        m_base_rates = ProcessDoubleVector(dict.at("base-rates"));

    if (dict.count("initial-taus") == 1)
        m_initial_taus = ProcessDoubleVector(dict.at("initial-taus"));

    // m_birth_times must be updated after m_g1_durations and m_g2_durations
    // are updated
    if (dict.count("birth-times") == 1)
        m_birth_times = ProcessDoubleVector(dict.at("birth-times"));
}

std::string DeathClockParameters::Print() const
{
    std::stringstream print_sstrm;

    print_sstrm << "output-directory = " << m_output_directory << std::endl;
    print_sstrm << "simulation-id = " << m_simulation_id << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "log-level = " << m_log_level << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "seed = " << m_seed << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "simulation-time = " << m_simulation_time << std::endl;
    print_sstrm << "dt = " << m_dt << std::endl;
    print_sstrm << "sampling-timestep-multiple = " <<
        m_sampling_timestep_multiple << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "cell-rearrangement-threshold = " <<
        m_cell_rearrangement_threshold << std::endl;
    print_sstrm << "cell-rearrangement-ratio = " << m_cell_rearrangement_ratio
        << std::endl;
    print_sstrm << "t2-threshold = " << m_t2_threshold << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "min-cell-count = " << m_min_cell_count << std::endl;
    print_sstrm << "max-cell-count = " << m_max_cell_count << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "min-cell-count-for-ancestor = " <<
        PrintUnsignedToUnsignedMap(m_min_cell_count_for_ancestor) << std::endl;
    print_sstrm << "max-cell-count-for-ancestor = " <<
        PrintUnsignedToUnsignedMap(m_max_cell_count_for_ancestor) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "num-cells-across = " << m_num_cells_across << std::endl;
    print_sstrm << "num-cells-up = " << m_num_cells_up << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "ancestors = " << PrintUnsignedVector(m_ancestors) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "force-model = " << m_force_model << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "boundary-tension-parameter = " << m_boundary_tension_parameter << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "cell-cycle-models = " << FormatUnsignedToCellCycleModel(PrintUnsignedVector(m_cell_cycle_models)) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "reference-target-areas = " << PrintDoubleVector(m_reference_target_areas) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "g1-durations = " << PrintDoubleVector(m_g1_durations) << std::endl;
    print_sstrm << "g2-durations = " << PrintDoubleVector(m_g2_durations) << std::endl;
    print_sstrm << "uniform-g1-ranges = " << PrintDoubleVector(m_uniform_g1_ranges) << std::endl;
    print_sstrm << "max-transit-generations = " << PrintUnsignedVector(m_max_transit_generations) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "apoptosis-times = " << PrintDoubleVector(m_apoptosis_times) << std::endl;
    print_sstrm << "death-thresholds = " << PrintDoubleVector(m_death_thresholds) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "base-rates = " << PrintDoubleVector(m_base_rates) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "neighbours-constants = " << PrintDoubleVector(m_neighbours_constants) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "neighbours-g2-constants = " << PrintDoubleVector(m_neighbours_g2_constants) << std::endl;
    print_sstrm << "neighbours-apoptosis-constants = " << PrintDoubleVector(m_neighbours_apoptosis_constants) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "normalised-neighbours-g2-constants = " << PrintDoubleVector(m_normalised_neighbours_g2_constants) << std::endl;
    print_sstrm << "normalised-neighbours-apoptosis-constants = " << PrintDoubleVector(m_normalised_neighbours_apoptosis_constants) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "shared-edge-constants = " << PrintDoubleVector(m_shared_edge_constants) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "shared-edge-g2-constants = " << PrintDoubleVector(m_shared_edge_g2_constants) << std::endl;
    print_sstrm << "shared-edge-apoptosis-constants = " << PrintDoubleVector(m_shared_edge_apoptosis_constants) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "normalised-shared-edge-g2-constants = " << PrintDoubleVector(m_normalised_shared_edge_g2_constants) << std::endl;
    print_sstrm << "normalised-shared-edge-apoptosis-constants = " << PrintDoubleVector(m_normalised_shared_edge_apoptosis_constants) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "initial-taus = " << PrintDoubleVector(m_initial_taus) << std::endl;

    print_sstrm << std::endl;
    print_sstrm << "birth-times = " << PrintDoubleVector(m_birth_times) << std::endl;

    return print_sstrm.str();
}

std::string DeathClockParameters::Help()
{
    std::string help_msg("Valid keys:\n");

    help_msg += "\n";
    help_msg += "  output-directory                         (default: death-clock)\n";
    help_msg += "  simulation-id                            (default: 0)\n";

    help_msg += "\n";
    help_msg += "  log-level                                (default: 1)\n";

    help_msg += "\n";
    help_msg += "  seed                                     (default: 0 (-> random))\n";

    help_msg += "\n";
    help_msg += "  simulation-time                          (default: 100)\n";
    help_msg += "  dt                                       (default: 0.01)\n";
    help_msg += "  sampling-timestep-multiple               (default: 100)\n";

    help_msg += "\n";
    help_msg += "  cell-rearrangement-threshold             (default: 0.01)\n";
    help_msg += "  cell-rearrangement-ratio                 (default: 1.5)\n";
    help_msg += "  t2-threshold                             (default: 0.001)\n";

    help_msg += "\n";
    help_msg += "  min-cell-count                           (default: 0)\n";
    help_msg += "  max-cell-count                           (default: inf)\n";

    help_msg += "\n";
    help_msg += "  min-cell-count-for-ancestor              (default: {})\n";
    help_msg += "  max-cell-count-for-ancestor              (default: {})\n";

    help_msg += "\n";
    help_msg += "  num-cells-across                         (default: 4)\n";
    help_msg += "  num-cells-up                             (default: 4)\n";

    help_msg += "\n";
    help_msg += "  ancestors                                (default: 0, 1, 2, ...)\n";

    help_msg += "\n";
    help_msg += "  force-model (farhadifar|nagai-honda)     (default: farhadifar)\n";

    help_msg += "\n";
    help_msg += "  boundary-tension-parameter               (default: 0.12)\n";

    help_msg += "\n";
    help_msg += "  cell-cycle-models (exponential|uniform)  (default: exponential)\n";

    help_msg += "\n";
    help_msg += "  reference-target-areas                   (default: 1.0)\n";

    help_msg += "\n";
    help_msg += "  g1-durations                             (default: 50.0)\n";
    help_msg += "  g2-durations                             (default: 50.0)\n";
    help_msg += "  uniform-g1-ranges                        (default: 50.0)\n";

    help_msg += "\n";
    help_msg += "  max-transit-generations                  (default: inf)\n";

    help_msg += "\n";
    help_msg += "  apoptosis-times                          (default: 0.25)\n";
    help_msg += "  death-thresholds                         (default: 100.0)\n";

    help_msg += "\n";
    help_msg += "  base-rates                               (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  neighbours-constants                     (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  neighbours-g2-constants                  (default: 0.0)\n";
    help_msg += "  neighbours-apoptosis-constants           (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  normalised-neighbours-g2-constants       (default: 6.0)\n";
    help_msg += "  normalised-neighbours-apoptosis-constants  (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  shared-edge-constants                     (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  shared-edge-g2-constants                  (default: 0.0)\n";
    help_msg += "  shared-edge-apoptosis-constants           (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  normalised-shared-edge-g2-constants       (default: 0.0)\n";
    help_msg += "  normalised-shared-edge-apoptosis-constants  (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  initial-taus                             (default: 0.0)\n";

    help_msg += "\n";
    help_msg += "  birth-times                              (default: random)\n";

    return help_msg;
}

std::vector<double> DeathClockParameters::ProcessDoubleVector(const std::string& input) const
{
    // Sanity check, m_initial_num_cells must be > 0 and
    // = m_num_cells_across * m_num_cells_up
    assert(m_initial_num_cells > 0);
    assert(m_initial_num_cells == m_num_cells_across * m_num_cells_up);

    // Read vector
    std::vector<double> vec = ReadDoubleVector(input);

    // Fill in based on number of elements in vector
    if (vec.size() == 1)
    {
        return std::vector<double>(m_initial_num_cells, vec[0]);
    }
    else if (vec.size() == m_initial_num_cells)
    {
        return vec;
    }
    else
    {
        std::string error_msg;
        error_msg += "vec.size() (";
        error_msg += std::to_string(vec.size());
        error_msg +=  ") is not equal to 1 or m_initial_num_cells (";
        error_msg += std::to_string(m_initial_num_cells);
        error_msg += ")";
        EXCEPTION(error_msg);
    }

    NEVER_REACHED;
}

std::vector<unsigned> DeathClockParameters::ProcessUnsignedVector(const std::string& input) const
{
    // Sanity check, m_initial_num_cells must be > 0 and
    // = m_num_cells_across * m_num_cells_up
    assert(m_initial_num_cells > 0);
    assert(m_initial_num_cells == m_num_cells_across * m_num_cells_up);

    // Read vector
    std::vector<unsigned> vec = ReadUnsignedVector(input);

    // Fill in based on number of elements in vector
    if (vec.size() == 1)
    {
        return std::vector<unsigned>(m_initial_num_cells, vec[0]);
    }
    else if (vec.size() == m_initial_num_cells)
    {
        return vec;
    }
    else
    {
        std::string error_msg;
        error_msg += "vec.size() (";
        error_msg += std::to_string(vec.size());
        error_msg += ") is not equal to 1 or m_initial_num_cells (";
        error_msg += std::to_string(m_initial_num_cells);
        error_msg += ")";
        EXCEPTION(error_msg);
    }

    NEVER_REACHED;
}

std::vector<double> DeathClockParameters::ReadDoubleVector(const std::string& input)
{
    std::vector<double> output;

    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ','))
    {
        output.push_back(std::stod(token));
    }

    return output;
}

std::vector<unsigned> DeathClockParameters::ReadUnsignedVector(const std::string& input)
{
    std::vector<unsigned> output;

    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ','))
    {
        output.push_back(std::stoul(token));
    }

    return output;
}

std::map<unsigned, unsigned> DeathClockParameters::ReadUnsignedToUnsignedMap(const std::string& input)
{
    std::map<unsigned, unsigned> output;

    if (input.front() != '{' || input.back() != '}')
    {
        std::string error_msg;
        error_msg += "Could not parse unsigned to unsigned map from string '";
        error_msg += input;
        error_msg += "'";
        error_msg += '\n';
        error_msg += "Correct format is '{KEY:VAL, KEY:VAL, ..., KEY:VAL}'";
        EXCEPTION(error_msg);
    }

    std::stringstream ss(input.substr(1, input.size() - 2));
    std::string token;

    while (std::getline(ss, token, ','))
    {
        auto separator_pos = token.find_first_of(':');

        if (separator_pos == std::string::npos)
        {
            std::string error_msg;
            error_msg += "Could not parse unsigned to unsigned map from string '";
            error_msg += input;
            error_msg += "'";
            error_msg += '\n';
            error_msg += "Correct format is '{KEY:VAL, KEY:VAL, ..., KEY:VAL}'";
            EXCEPTION(error_msg);
        }

        auto key = std::stoul(token.substr(0, separator_pos));
        auto val = std::stoul(token.substr(separator_pos + 1, std::string::npos));

        output[key] = val;
    }

    return output;
}

std::string DeathClockParameters::PrintDoubleVector(const std::vector<double>& input)
{
    std::string output;
    for (double num : input)
    {
        output += std::to_string(num);
        output += ", ";
    }

    // Delete last comma and space
    if (input.size() > 0)
    {
        output.pop_back();
        output.pop_back();
    }

    return output;
}

std::string DeathClockParameters::PrintUnsignedVector(const std::vector<unsigned>& input)
{
    std::string output;
    for (unsigned num : input)
    {
        output += std::to_string(num);
        output += ", ";
    }

    // Delete last comma and space
    if (input.size() > 0)
    {
        output.pop_back();
        output.pop_back();
    }

    return output;
}

std::string DeathClockParameters::PrintUnsignedToUnsignedMap(const
        std::map<unsigned, unsigned>& input)
{
    std::string output;

    output += '{';

    for (const auto& item : input)
    {
        output += std::to_string(item.first);
        output += ':';
        output += std::to_string(item.second);
        output += ", ";
    }

    // Delete last comma and space
    if (input.size() > 0)
    {
        output.pop_back();
        output.pop_back();
    }

    output += '}';

    return output;
}

std::string DeathClockParameters::FormatUnsignedToCellCycleModel(const std::string& input)
{
    std::string output;

    std::stringstream ss(input);
    std::string token;

    bool at_least_one_token = false;
    while (std::getline(ss, token, ','))
    {
        at_least_one_token = true;
        CellCycleModel ccm = static_cast<CellCycleModel>(std::stoul(token));

        output += CellCycleModelToString(ccm);
        output += ", ";
    }

    //Delete last comma and space
    if (at_least_one_token)
    {
        output.pop_back();
        output.pop_back();
    }

    return output;
}

std::string DeathClockParameters::FormatCellCycleModelToUnsigned(const std::string& input)
{
    std::string output;

    std::stringstream ss(input);
    std::string token;

    bool at_least_one_token = false;
    while (std::getline(ss, token, ','))
    {
        at_least_one_token = true;
        remove_leading_trailing_whitespace(token);
        CellCycleModel ccm = StringToCellCycleModel(token);

        output += std::to_string(static_cast<unsigned>(ccm));
        output += ", ";
    }

    //Delete last comma and space
    if (at_least_one_token)
    {
        output.pop_back();
        output.pop_back();
    }

    return output;
}

bool DeathClockParameters::IsValid() const
{
    // Criteria for IsValid:
    // - m_initial_num_cells = m_num_cells_across * m_num_cells_up
    // - all vector quantities have size m_initial_num_cells
    // - m_birth_times are in the range [-(m_g1_durations + m_g2_durations), 0]

    if (m_initial_num_cells != m_num_cells_across * m_num_cells_up)
    {
        return false;
    }

    if (   m_ancestors.size() != m_initial_num_cells
        || m_cell_cycle_models.size() != m_initial_num_cells
        || m_reference_target_areas.size() != m_initial_num_cells
        || m_g1_durations.size() != m_initial_num_cells
        || m_g2_durations.size() != m_initial_num_cells
        || m_uniform_g1_ranges.size() != m_initial_num_cells
        || m_max_transit_generations.size() != m_initial_num_cells
        || m_apoptosis_times.size() != m_initial_num_cells
        || m_death_thresholds.size() != m_initial_num_cells
        || m_neighbours_constants.size() != m_initial_num_cells
        || m_neighbours_g2_constants.size() != m_initial_num_cells
        || m_neighbours_apoptosis_constants.size() != m_initial_num_cells
        || m_normalised_neighbours_g2_constants.size() != m_initial_num_cells
        || m_normalised_neighbours_apoptosis_constants.size() != m_initial_num_cells
        || m_shared_edge_constants.size() != m_initial_num_cells
        || m_shared_edge_g2_constants.size() != m_initial_num_cells
        || m_shared_edge_apoptosis_constants.size() != m_initial_num_cells
        || m_normalised_shared_edge_g2_constants.size() != m_initial_num_cells
        || m_normalised_shared_edge_apoptosis_constants.size() != m_initial_num_cells
        || m_base_rates.size() != m_initial_num_cells
        || m_initial_taus.size() != m_initial_num_cells
        || m_birth_times.size() != m_initial_num_cells)
    {
        return false;
    }

    for (unsigned i = 0; i < m_initial_num_cells; i++)
    {
        double birth_time = m_birth_times[i];
        double min_age = - (m_g1_durations[i] + m_g2_durations[i]);

        if (0 < birth_time || birth_time < min_age)
        {
            return false;
        }
    }

    return true;
}
