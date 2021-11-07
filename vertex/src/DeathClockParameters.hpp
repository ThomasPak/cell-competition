#ifndef DEATHCLOCKPARAMETERS_HPP
#define DEATHCLOCKPARAMETERS_HPP

#include <string>
#include <vector>
#include <map>

#include "Dictionary.hpp"

enum CellCycleModel : unsigned { exponential, uniform };

CellCycleModel StringToCellCycleModel(const std::string& input);
std::string CellCycleModelToString(CellCycleModel ccm);

/** Struct that contains parameters for DeathClockSimulation().
 *
 */
struct DeathClockParameters
{

    // Default constructor
    DeathClockParameters();

    /** Output directory.  Defaults to "death-clock". */
    std::string m_output_directory = "death-clock";
    std::string m_simulation_id = "0";

    // Log level
    unsigned m_log_level = 1;

    // Define simulation parameters
    unsigned m_seed = 0;

    double m_simulation_time = 100;
    double m_dt = 0.01;
    unsigned m_sampling_timestep_multiple = 100;

    double m_cell_rearrangement_threshold = 0.01;
    double m_cell_rearrangement_ratio = 1.5;
    double m_t2_threshold = 0.001;

    unsigned m_min_cell_count = 0u;
    unsigned m_max_cell_count = static_cast<unsigned>(-1);

    std::map<unsigned, unsigned> m_min_cell_count_for_ancestor;
    std::map<unsigned, unsigned> m_max_cell_count_for_ancestor;

    unsigned m_num_cells_across = 4;
    unsigned m_num_cells_up = 4;

    unsigned m_initial_num_cells = m_num_cells_across * m_num_cells_up;

    // Ancestors will be initialised in the default constructor
    std::vector<unsigned> m_ancestors;

    // Define model parameters
    std::string m_force_model = "farhadifar";

    // Farhadifar force parameters
    double m_boundary_tension_parameter = 0.12;

    std::vector<unsigned> m_cell_cycle_models = std::vector<unsigned>(m_initial_num_cells, exponential);

    std::vector<double> m_reference_target_areas = std::vector<double>(m_initial_num_cells, 1.0);

    std::vector<double> m_g1_durations = std::vector<double>(m_initial_num_cells, 50.0);
    std::vector<double> m_g2_durations = std::vector<double>(m_initial_num_cells, 50.0);
    std::vector<double> m_uniform_g1_ranges = std::vector<double>(m_initial_num_cells, 10.0);
    std::vector<unsigned> m_max_transit_generations = std::vector<unsigned>(m_initial_num_cells, static_cast<unsigned>(-1));

    std::vector<double> m_apoptosis_times = std::vector<double>(m_initial_num_cells, 0.25);
    std::vector<double> m_death_thresholds = std::vector<double>(m_initial_num_cells, 100.0);

    std::vector<double> m_neighbours_constants = std::vector<double>(m_initial_num_cells, 0.0);

    std::vector<double> m_neighbours_g2_constants = std::vector<double>(m_initial_num_cells, 0.0);
    std::vector<double> m_neighbours_apoptosis_constants = std::vector<double>(m_initial_num_cells, 0.0);

    std::vector<double> m_normalised_neighbours_g2_constants = std::vector<double>(m_initial_num_cells, 6.0);
    std::vector<double> m_normalised_neighbours_apoptosis_constants = std::vector<double>(m_initial_num_cells, 0.0);

    std::vector<double> m_shared_edge_constants = std::vector<double>(m_initial_num_cells, 0.0);

    std::vector<double> m_shared_edge_g2_constants = std::vector<double>(m_initial_num_cells, 0.0);
    std::vector<double> m_shared_edge_apoptosis_constants = std::vector<double>(m_initial_num_cells, 0.0);

    std::vector<double> m_normalised_shared_edge_g2_constants = std::vector<double>(m_initial_num_cells, 0.0);
    std::vector<double> m_normalised_shared_edge_apoptosis_constants = std::vector<double>(m_initial_num_cells, 0.0);

    std::vector<double> m_base_rates = std::vector<double>(m_initial_num_cells, 0.0);

    // Define initial state
    std::vector<double> m_initial_taus = std::vector<double>(m_initial_num_cells, 0.0);

    // Birth times will be initialised in the default constructor
    std::vector<double> m_birth_times;

    /** Update values given by a Dictionary.
     *
     * @param dict  Dictionary containing updated values.
     */
    void Update(const Dictionary& dict);

    /** Check whether struct is valid */
    bool IsValid() const;

    /** Print parameter values.
     *
     * @return string containing printed parameter values.
     */
    std::string Print() const;

    /** Process vector of doubles
     *
     * @return vector of double values.
     */
    std::vector<double> ProcessDoubleVector(const std::string& input) const;

    /** Read vector of unsigneds
     *
     * @return vector of unsigned values.
     */
    std::vector<unsigned> ProcessUnsignedVector(const std::string& input) const;

    /** Generate Help message
     *
     * @return help string.
     */
    static std::string Help();

    /** Read vector of doubles
     *
     * @return vector of double values.
     */
    static std::vector<double> ReadDoubleVector(const std::string& input);

    /** Read vector of unsigneds
     *
     * @return vector of unsigned values.
     */
    static std::vector<unsigned> ReadUnsignedVector(const std::string& input);

    /** Read map of unsigned to unsigned values
     *
     * @return map of unsigned to unsigned values.
     */
    static std::map<unsigned, unsigned> ReadUnsignedToUnsignedMap(const std::string& input);

    /** Print vector of doubles
     *
     * @return string representing vector of double values
     */
    static std::string PrintDoubleVector(const std::vector<double>& input);

    /** Print vector of unsigneds
     *
     * @return string representing vector of double values
     */
    static std::string PrintUnsignedVector(const std::vector<unsigned>& input);

    /** Print map of unsigned to unsigned values
     *
     * @return string representing map of unsigned to unsigned values.
     */
    static std::string PrintUnsignedToUnsignedMap(const std::map<unsigned, unsigned>& input);

    /** Print vector of cell cycle models
     *
     * @return string representing vector of cell cycle models
     */
    static std::string PrintCellCycleModelVector(const std::vector<unsigned>& input);

    /** Convert unsigned string to cell cycle model string
     */
    static std::string FormatUnsignedToCellCycleModel(const std::string& input);

    /** Convert cell cycle model string to unsigned string
     */
    static std::string FormatCellCycleModelToUnsigned(const std::string& input);

};

#endif // DEATHCLOCKPARAMETERS_HPP
