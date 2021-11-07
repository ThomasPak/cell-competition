#ifndef TESTDEATHCLOCKODESYSTEM_HPP_
#define TESTDEATHCLOCKODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "DeathClockOdeSystem.hpp"
#include "CvodeAdaptor.hpp"
#include "Timer.hpp"
#include "OutputFileHandler.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestDeathClockOdeSystem : public CxxTest::TestSuite
{
public:

    void TestDeathClockOdeSystemSetup()
    {
#ifdef CHASTE_CVODE
        DeathClockOdeSystem ode_system;

        // Set parameters
        ode_system.SetParameter("number of neighbours", 5);
        ode_system.SetParameter("number of neighbours g2", 2);
        ode_system.SetParameter("number of neighbours apoptosis", 1);

        ode_system.SetParameter("shared edge length", 50);
        ode_system.SetParameter("shared edge length g2", 20);
        ode_system.SetParameter("shared edge length apoptosis", 10);

        ode_system.SetParameter("base rate", 0.75);

        ode_system.SetParameter("neighbours constant", .76);
        ode_system.SetParameter("neighbours g2 constant", .77);
        ode_system.SetParameter("neighbours apoptosis constant", .78);
        ode_system.SetParameter("normalised neighbours g2 constant", .79);
        ode_system.SetParameter("normalised neighbours apoptosis constant", .80);

        ode_system.SetParameter("shared edge constant", .076);
        ode_system.SetParameter("shared edge g2 constant", .077);
        ode_system.SetParameter("shared edge apoptosis constant", .078);
        ode_system.SetParameter("normalised shared edge g2 constant", .079);
        ode_system.SetParameter("normalised shared edge apoptosis constant", .080);

        double h_value = 1.0;
        CvodeAdaptor cvode_solver;
        OdeSolution solutions;

        // Check initial conditions
        std::vector<double> initial_conditions = ode_system.GetInitialConditions();

        TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-6);

        Timer::Reset();
        solutions = cvode_solver.Solve(&ode_system, initial_conditions, 0.0, 100.0, h_value, h_value);
        Timer::Print("1. Cvode for 100 hours");

        // Test solutions are OK for a small time increase...
        int end = solutions.rGetSolutions().size() - 1;
        // Tests the simulation is ending at the right time
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 100, 1e-2);

        // Decent results
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 1351.36, 1e-4);
#else
        std::cout << "CVODE is not enabled. " << std::endl;
        std::cout << "If required please install and alter your hostconfig "
            "settings to switch on chaste support." << std::endl;
#endif //CHASTE_CVODE
    }

    void TestArchiving()
    {
#ifdef CHASTE_CVODE
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "death_clock_ode.arch";

        {
            std::vector<double> state_variables;
            state_variables.push_back(3.0);

            DeathClockOdeSystem ode_system(state_variables);

            ode_system.SetDefaultInitialCondition(0, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();

            // These are the initial conditions hard-coded in the constructor.
            TS_ASSERT_EQUALS(initial_conditions.size(), 1u);
            TS_ASSERT_DELTA(initial_conditions[0], 3.25, 1e-6);

            ode_system.SetParameter("base rate", 10.0);

            ode_system.SetParameter("neighbours constant", 10.1);
            ode_system.SetParameter("neighbours g2 constant", 10.2);
            ode_system.SetParameter("neighbours apoptosis constant", 10.3);
            ode_system.SetParameter("normalised neighbours g2 constant", 10.4);
            ode_system.SetParameter("normalised neighbours apoptosis constant", 10.5);

            ode_system.SetParameter("shared edge constant", 20.1);
            ode_system.SetParameter("shared edge g2 constant", 20.2);
            ode_system.SetParameter("shared edge apoptosis constant", 20.3);
            ode_system.SetParameter("normalised shared edge g2 constant", 20.4);
            ode_system.SetParameter("normalised shared edge apoptosis constant", 20.5);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();

            TS_ASSERT_EQUALS(initial_conditions.size(), 1u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.0, 1e-6);

            double base_rate = p_ode_system->GetParameter("base rate");

            double neighbours_constant = p_ode_system->GetParameter("neighbours constant");
            double neighbours_g2_constant =
                p_ode_system->GetParameter("neighbours g2 constant");
            double neighbours_apoptosis_constant =
                p_ode_system->GetParameter("neighbours apoptosis constant");
            double normalised_neighbours_g2_constant =
                p_ode_system->GetParameter("normalised neighbours g2 constant");
            double normalised_neighbours_apoptosis_constant =
                p_ode_system->GetParameter("normalised neighbours apoptosis constant");

            double shared_edge_constant = p_ode_system->GetParameter("shared edge constant");
            double shared_edge_g2_constant =
                p_ode_system->GetParameter("shared edge g2 constant");
            double shared_edge_apoptosis_constant =
                p_ode_system->GetParameter("shared edge apoptosis constant");
            double normalised_shared_edge_g2_constant =
                p_ode_system->GetParameter("normalised shared edge g2 constant");
            double normalised_shared_edge_apoptosis_constant =
                p_ode_system->GetParameter("normalised shared edge apoptosis constant");

            TS_ASSERT_DELTA(base_rate, 10.0, 1e-3);

            TS_ASSERT_DELTA(neighbours_constant, 10.1, 1e-3);
            TS_ASSERT_DELTA(neighbours_g2_constant, 10.2, 1e-3);
            TS_ASSERT_DELTA(neighbours_apoptosis_constant, 10.3, 1e-3);
            TS_ASSERT_DELTA(normalised_neighbours_g2_constant, 10.4, 1e-3);
            TS_ASSERT_DELTA(normalised_neighbours_apoptosis_constant, 10.5, 1e-3);

            TS_ASSERT_DELTA(shared_edge_constant, 20.1, 1e-3);
            TS_ASSERT_DELTA(shared_edge_g2_constant, 20.2, 1e-3);
            TS_ASSERT_DELTA(shared_edge_apoptosis_constant, 20.3, 1e-3);
            TS_ASSERT_DELTA(normalised_shared_edge_g2_constant, 20.4, 1e-3);
            TS_ASSERT_DELTA(normalised_shared_edge_apoptosis_constant, 20.5, 1e-3);

            // Get state variables
            double var1 = p_ode_system->GetStateVariable(0);

            TS_ASSERT_DELTA(var1, 3.0, 1e-3);

            // Tidy up
            delete p_ode_system;
        }
#endif //CHASTE_CVODE
    }

    void TestSetStateVariables()
    {
#ifdef CHASTE_CVODE

        std::vector<double> state_vars;
        state_vars.push_back(1.0);
        DeathClockOdeSystem ode_system(state_vars);

        TS_ASSERT_EQUALS(ode_system.GetStateVariable(0), 1.0);

#endif //CHASTE_CVODE
   }
};

#endif // TESTDEATHCLOCKODESYSTEM_HPP_
