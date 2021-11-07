#ifndef TESTCELLCOMPETITION_HPP_
#define TESTCELLCOMPETITION_HPP_

#include <iostream>
#include <string>
#include <stdlib.h>
#include <string.h>

#include <cxxtest/TestSuite.h>

#include "ExecutableSupport.hpp"

#include "Dictionary.hpp"
#include "SingletonSupport.hpp"
#include "CellCompetition.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCellCompetition : public CxxTest::TestSuite
{
    public:

        void TestRunCellCompetitionSimulation()
        {
            // Setup parameters
            CellCompetitionParameters param;

            param.m_simulation_time = 1;
            param.m_simulation_time = 0.01;

            param.m_num_cells_across = 2;
            param.m_num_cells_up = 2;

            // Setup singletons
            SingletonSupport::SetupSingletons(1);

            // Run simulation
            RunCellCompetitionSimulation(param);

            // Destroy singletons
            SingletonSupport::DestroySingletons();
        }

        void TestCellCompetitionSimulatorHelp1()
        {
            // Setup arguments
            int argc = 2;
            char* argv[2];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("--help") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "--help");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    CellCompetitionSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestCellCompetitionSimulatorHelp2()
        {
            // Setup arguments
            int argc = 2;
            char* argv[2];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("-h") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "-h");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    CellCompetitionSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestCellCompetitionSimulatorInvalidArgument()
        {
            // Setup arguments
            int argc = 2;
            char* argv[2];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("blabla=") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "blabla=");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    CellCompetitionSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_ERROR)

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestCellCompetitionSimulatorRunWithArguments()
        {
            // Setup arguments
            int argc = 5;
            char* argv[5];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            argv[1] = static_cast<char*>(malloc(sizeof(char) * (strlen("simulation-time=1") + 1)));
            argv[2] = static_cast<char*>(malloc(sizeof(char) * (strlen("num-cells-across=2") + 1)));
            argv[3] = static_cast<char*>(malloc(sizeof(char) * (strlen("num-cells-up=2") + 1)));
            argv[4] = static_cast<char*>(malloc(sizeof(char) * (strlen("seed=1") + 1)));
            strcpy(argv[0], "test");
            strcpy(argv[1], "simulation-time=1");
            strcpy(argv[2], "num-cells-across=2");
            strcpy(argv[3], "num-cells-up=2");
            strcpy(argv[4], "seed=1");

            // Setup input and output string
            std::string input_string, output_string;

            // Run simulator
            TS_ASSERT_EQUALS(
                    CellCompetitionSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }

        void TestCellCompetitionSimulatorRunWithInput()
        {
            // Setup arguments
            int argc = 1;
            char* argv[1];
            argv[0] = static_cast<char*>(malloc(sizeof(char) * (strlen("test") + 1)));
            strcpy(argv[0], "test");

            // Setup input and output string
            std::string input_string, output_string;
            input_string += "simulation-time=1\n";
            input_string += "num-cells-across=2\n";
            input_string += "num-cells-up=2\n";
            input_string += "seed=1\n";

            // Run simulator
            TS_ASSERT_EQUALS(
                    CellCompetitionSimulator(argc, argv, input_string, output_string),
                    ExecutableSupport::EXIT_OK);

            TS_ASSERT_EQUALS(output_string, "");
        }
};

#endif // TESTCELLCOMPETITION_HPP_
