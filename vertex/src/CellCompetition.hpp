#ifndef CELLCOMPETITION_HPP_
#define CELLCOMPETITION_HPP_

#include <string>

#include "Dictionary.hpp"
#include "CellCompetitionParameters.hpp"

int CellCompetitionSimulator(int argc, char *argv[],
        const std::string& input_string, std::string& output_string);

void RunCellCompetitionSimulation(const CellCompetitionParameters& param);

#endif // CELLCOMPETITION_HPP_
