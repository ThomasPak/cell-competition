#include <iostream>
#include <string>

#include <unistd.h>

#include "ExecutableSupport.hpp"

#include "DeathClockParameters.hpp"
#include "DeathClockSimulator.hpp"

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
	ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    // Process arguments
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help" ))
    {
        std::cerr << "Usage: " << argv[0] << " [KEY=VALUE]...\n";
        std::cerr << std::endl;
        std::cerr << DeathClockParameters::Help() << std::endl;

        ExecutableSupport::FinalizePetsc();
        return ExecutableSupport::EXIT_OK;
    }

    // Prepare input and output strings
    std::string input_string, output_string;

    // Prompt for input
    if (isatty(STDIN_FILENO))
        std::cerr << "Enter parameters as (multiple non-empty lines allowed):\n"
            "KEY_1=VALUE_1\nKEY_2=VALUE_2\n...\nKEY_M=VALUE_M\n";

    std::string line;
    while (std::getline(std::cin, line))
    {
        input_string += line;
        input_string += '\n';
    }

    // Run simulation
    int exit_code = DeathClockSimulator(argc, argv, input_string,
            output_string);

    // Write output to stdout
    std::cout << output_string;

    // End by finalizing PETSc, and returning a suitable exit code.
    // 0 means 'no error'
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
