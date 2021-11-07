#ifndef TESTDEATHCLOCKSRNMODEL_HPP_
#define TESTDEATHCLOCKSRNMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "DeathClockSrnModel.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestDeathClockSrnModel : public AbstractCellBasedTestSuite
{
    public:

        void TestDeathClockSrnBaseRate()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set base rate to 3.0
            p_srn_model->SetBaseRate(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetBaseRate(), 3.0, 1e-4);

            // Set number of neighbours to 0.0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 30.5, 1e-4);
        }

        void TestDeathClockSrnNeighbours()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set neighbours constant to 3.0
            p_srn_model->SetNeighboursConstant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetNeighboursConstant(), 3.0, 1e-4);

            // Set number of neighbours artifically to 2.0
            p_cell->GetCellData()->SetItem("number of neighbours", 2);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 60.5, 1e-4);
        }

        void TestDeathClockSrnNeighboursG2()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set neighbours g2 constant to 3.0
            p_srn_model->SetNeighboursG2Constant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetNeighboursG2Constant(), 3.0, 1e-4);

            // Set number of neighbours G2 artifically to 2.0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 2.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 60.5, 1e-4);
        }

        void TestDeathClockSrnNeighboursApoptosis()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set neighbours apoptosis constant to 3.0
            p_srn_model->SetNeighboursApoptosisConstant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetNeighboursApoptosisConstant(), 3.0, 1e-4);

            // Set number of neighbours apoptosis artifically to 2.0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 2.0);

            // Set shared edge length to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 60.5, 1e-4);
        }

        void TestDeathClockSrnNormalisedNeighboursG2()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set normalised neighbours g2 constant to 3.0
            p_srn_model->SetNormalisedNeighboursG2Constant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetNormalisedNeighboursG2Constant(), 3.0, 1e-4);

            // Set number of neighbours G2 artifically to 2.0
            p_cell->GetCellData()->SetItem("number of neighbours", 5.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 2.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 12.5, 1e-4);
        }

        void TestDeathClockSrnNormalisedNeighboursApoptosis()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set normalised neighbours apoptosis constant to 3.0
            p_srn_model->SetNormalisedNeighboursApoptosisConstant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetNormalisedNeighboursApoptosisConstant(), 3.0, 1e-4);

            // Set number of neighbours apoptosis artifically to 2.0
            p_cell->GetCellData()->SetItem("number of neighbours", 5.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 2.0);

            // Set shared edge length to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 12.5, 1e-4);
        }

        void TestDeathClockSrnSharedEdge()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set shared edge constant to 3.0
            p_srn_model->SetSharedEdgeConstant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetSharedEdgeConstant(), 3.0, 1e-4);

            // Set number of neighbours to 0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length to artificially to 2.0
            p_cell->GetCellData()->SetItem("shared edge length", 2.5);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 75.5, 1e-4);
        }

        void TestDeathClockSrnSharedEdgeG2()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set shared edge g2 constant to 3.0
            p_srn_model->SetSharedEdgeG2Constant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetSharedEdgeG2Constant(), 3.0, 1e-4);

            // Set number of neighbours G2 to 0.0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length G2 to 2.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 2.5);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 75.5, 1e-4);
        }

        void TestDeathClockSrnSharedEdgeApoptosis()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set shared edge apoptosis constant to 3.0
            p_srn_model->SetSharedEdgeApoptosisConstant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetSharedEdgeApoptosisConstant(), 3.0, 1e-4);

            // Set number of neighbours apoptosis to 0.0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length apoptosis to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 2.5);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 75.5, 1e-4);
        }

        void TestDeathClockSrnNormalisedSharedEdgeG2()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set normalised shared edge g2 constant to 3.0
            p_srn_model->SetNormalisedSharedEdgeG2Constant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetNormalisedSharedEdgeG2Constant(), 3.0, 1e-4);

            // Set number of neighbours G2 to 0.0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length G2 to 2.0
            p_cell->GetCellData()->SetItem("shared edge length", 5.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 2.5);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 0.0);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 15.5, 1e-4);
        }

        void TestDeathClockSrnNormalisedSharedEdgeApoptosis()
        {
            TS_ASSERT_THROWS_NOTHING(DeathClockSrnModel srn_model);

            DeathClockSrnModel* p_srn_model = new DeathClockSrnModel();

            // Create a vector of initial conditions
            std::vector<double> starter_conditions;
            starter_conditions.push_back(0.5);
            p_srn_model->SetInitialConditions(starter_conditions);

            UniformG1GenerationalCellCycleModel* p_cc_model = new
                UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model,
                        false, CellPropertyCollection()));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            // Now updated to initial conditions
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 0.5, 1e-4);

            // Set normalised shared edge apoptosis constant to 3.0
            p_srn_model->SetNormalisedSharedEdgeApoptosisConstant(3.0);
            TS_ASSERT_DELTA(p_srn_model->GetNormalisedSharedEdgeApoptosisConstant(), 3.0, 1e-4);

            // Set number of neighbours apoptosis to 0.0
            p_cell->GetCellData()->SetItem("number of neighbours", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours g2", 0.0);
            p_cell->GetCellData()->SetItem("number of neighbours apoptosis", 0.0);

            // Set shared edge length to 0.0
            p_cell->GetCellData()->SetItem("shared edge length", 5.0);
            p_cell->GetCellData()->SetItem("shared edge length g2", 0.0);
            p_cell->GetCellData()->SetItem("shared edge length apoptosis", 2.5);

            // Now update the SRN
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            unsigned num_steps = 100;
            double end_time = 10.0;
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            // Test evolved as expected
            TS_ASSERT_DELTA(p_srn_model->GetTau(), 15.5, 1e-4);
        }


        void TestDeathClockSrnCreateCopy()
        {
            // Test with DeathClockSrnModel
            DeathClockSrnModel* p_model= new DeathClockSrnModel;

            // Set ODE system
            std::vector<double> state_variables;
            state_variables.push_back(2.0);
            p_model->SetOdeSystem(new DeathClockOdeSystem(state_variables));

            p_model->SetInitialConditions(state_variables);

            // Create a copy
            DeathClockSrnModel* p_model2 = static_cast<DeathClockSrnModel*>
                (p_model->CreateSrnModel());

            // Check correct initializations
            TS_ASSERT_EQUALS(p_model2->GetTau(), 2.0);

            // Destroy models
            delete p_model;
            delete p_model2;
        }

        void TestArchiveDeathClockSrnModel()
        {
            OutputFileHandler handler("archive", false);
            std::string archive_filename = handler.GetOutputDirectoryFullPath()
                + "death_clock_srn.arch";

            // Create an output archive
            {
                SimulationTime* p_simulation_time = SimulationTime::Instance();
                p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

                UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

                // As usual, we archive via a pointer to the most abstract class possible
                AbstractSrnModel* p_srn_model = new DeathClockSrnModel;

                MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
                MAKE_PTR(TransitCellProliferativeType, p_transit_type);

                // We must create a cell to be able to initialise the cell srn model's ODE system
                CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->InitialiseCellCycleModel();
                p_cell->InitialiseSrnModel();
                p_cell->SetBirthTime(0.0);

                static_cast<DeathClockSrnModel*>(p_srn_model)->SetBaseRate(0.5);

                static_cast<DeathClockSrnModel*>(p_srn_model)->SetNeighboursConstant(0.6);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetNeighboursG2Constant(0.7);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetNeighboursApoptosisConstant(0.8);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetNormalisedNeighboursG2Constant(0.9);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetNormalisedNeighboursApoptosisConstant(1.0);

                static_cast<DeathClockSrnModel*>(p_srn_model)->SetSharedEdgeConstant(1.6);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetSharedEdgeG2Constant(1.7);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetSharedEdgeApoptosisConstant(1.8);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetNormalisedSharedEdgeG2Constant(1.9);
                static_cast<DeathClockSrnModel*>(p_srn_model)->SetNormalisedSharedEdgeApoptosisConstant(2.0);

                std::ofstream ofs(archive_filename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                output_arch << p_srn_model;

                // Note that here, deletion of the cell-cycle model and srn is handled by the cell destructor
                SimulationTime::Destroy();
            }

            {
                // We must set SimulationTime::mStartTime here to avoid tripping an assertion
                SimulationTime::Instance()->SetStartTime(0.0);

                AbstractSrnModel* p_srn_model;

                std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
                boost::archive::text_iarchive input_arch(ifs);

                input_arch >> p_srn_model;

                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetBaseRate(), 0.5, 1e-12);

                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetNeighboursConstant(), 0.6, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetNeighboursG2Constant(), 0.7, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetNeighboursApoptosisConstant(), 0.8, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetNormalisedNeighboursG2Constant(), 0.9, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetNormalisedNeighboursApoptosisConstant(), 1.0, 1e-12);

                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetSharedEdgeConstant(), 1.6, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetSharedEdgeG2Constant(), 1.7, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetSharedEdgeApoptosisConstant(), 1.8, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetNormalisedSharedEdgeG2Constant(), 1.9, 1e-12);
                TS_ASSERT_DELTA(static_cast<DeathClockSrnModel*>(p_srn_model)->GetNormalisedSharedEdgeApoptosisConstant(), 2.0, 1e-12);

                delete p_srn_model;
            }
        }
};

#endif // TESTDEATHCLOCKSRNMODEL_HPP_
