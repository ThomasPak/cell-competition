#include <iostream>
#include <string>
#include <sstream>

#include "ExecutableSupport.hpp"
#include "Exception.hpp"

#include "SingletonSupport.hpp"
#include "PakmanMPIWorker.hpp"
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

    // Create MPI Worker
    PakmanMPIWorker worker(&DeathClockSimulator);

    // Run MPI Worker
    worker.run(argc, argv);

    // End by finalizing PETSc, and returning a suitable exit code.
    // 0 means 'no error'
    ExecutableSupport::FinalizePetsc();
    return ExecutableSupport::EXIT_OK;
}
