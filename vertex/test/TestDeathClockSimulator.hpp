#ifndef TESTDEATHCLOCKSIMULATOR_HPP
#define TESTDEATHCLOCKSIMULATOR_HPP

#include <cxxtest/TestSuite.h>

#include "ExecutableSupport.hpp"

#include "Dictionary.hpp"
#include "SingletonSupport.hpp"
#include "DeathClockSimulator.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDeathClock : public CxxTest::TestSuite
{
    public:

        void TestRunDeathClockSimulation()
        {
            // Setup parameters
            DeathClockParameters param;

            Dictionary dict  = create_dictionary(""
            "simulation-id = 0\n"
            "simulation-time = 1\n"
            "num-cells-across = 2\n"
            "num-cells-up = 2\n"
            );

            param.Update(dict);

            // Setup singletons
            SingletonSupport::SetupSingletons(1);

            // Run simulation
            RunDeathClockSimulation(param);

            // Destroy singletons
            SingletonSupport::DestroySingletons();
        }

        void TestDeathClockSimulatorHelp1()
        {
            // Setup arguments
            constexpr int argc = 2;
            char* argv[argc];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("--help") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "--help");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    DeathClockSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestDeathClockSimulatorHelp2()
        {
            // Setup arguments
            constexpr int argc = 2;
            char* argv[argc];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("-h") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "-h");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    DeathClockSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestDeathClockSimulatorInvalidArgument()
        {
            // Setup arguments
            constexpr int argc = 2;
            char* argv[argc];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("blabla=") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "blabla=");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    DeathClockSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_ERROR)

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestDeathClockSimulatorRunWithArguments()
        {
            // Setup arguments
            constexpr int argc = 6;
            char* argv[argc];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("simulation-id=1") + 1)));
            argv[2] = static_cast<char*>(malloc(sizeof(char) * (strlen("simulation-time=1") + 1)));
            argv[3] = static_cast<char*>(malloc(sizeof(char) * (strlen("num-cells-across=3") + 1)));
            argv[4] = static_cast<char*>(malloc(sizeof(char) * (strlen("num-cells-up=2") + 1)));
            argv[5] = static_cast<char*>(malloc(sizeof(char) * (strlen("seed=1") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "simulation-id=1");
            strcpy(argv[2], "simulation-time=1");
            strcpy(argv[3], "num-cells-across=3");
            strcpy(argv[4], "num-cells-up=2");
            strcpy(argv[5], "seed=1");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    DeathClockSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestDeathClockSimulatorRunWithInput()
        {
            // Setup arguments
            int argc = 1;
            char* argv[1];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            strcpy(argv[0], "test");

            // Setup input and output string
            std::string input_string, output_string;
            input_string += "simulation-id=2\n";
            input_string += "simulation-time=1\n";
            input_string += "num-cells-across=4\n";
            input_string += "num-cells-up=2\n";
            input_string += "seed=1\n";

            // Run simulator
            TS_ASSERT_EQUALS(
                    DeathClockSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }
};

#endif // TESTDEATHCLOCKSIMULATOR_HPP
