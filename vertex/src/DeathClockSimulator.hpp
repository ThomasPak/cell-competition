#ifndef DEATHCLOCKSIMULATOR_HPP
#define DEATHCLOCKSIMULATOR_HPP

#include <string>

#include "Dictionary.hpp"
#include "DeathClockParameters.hpp"

int DeathClockSimulator(int argc, char *argv[],
        const std::string& input_string, std::string& output_string);

void RunDeathClockSimulation(const DeathClockParameters& param);

#endif // DEATHCLOCKSIMULATOR_HPP
