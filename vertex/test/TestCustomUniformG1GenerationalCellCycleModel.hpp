#ifndef TESTCUSTOMUNIFORMG1GENERATIONALCELLCYCLEMODEL_HPP
#define TESTCUSTOMUNIFORMG1GENERATIONALCELLCYCLEMODEL_HPP

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FileComparison.hpp"
#include "OutputFileHandler.hpp"
#include "CustomUniformG1GenerationalCellCycleModel.hpp"

#include "RandomNumberGenerator.hpp"

class TestCustomUniformG1GenerationalCellCycleModel : public AbstractCellBasedTestSuite
{
public:

    void TestCustomUniformG1GenerationalCellCycleModel1()
    {
        TS_ASSERT_THROWS_NOTHING(CustomUniformG1GenerationalCellCycleModel cell_model3);

        CustomUniformG1GenerationalCellCycleModel* p_stem_model = new CustomUniformG1GenerationalCellCycleModel;


        // Change G1 duration and range for this model
        p_stem_model->SetStemCellG1Duration(3.0);
        p_stem_model->SetRange(4.0);

        CustomUniformG1GenerationalCellCycleModel* p_transit_model = new CustomUniformG1GenerationalCellCycleModel;

        // Change G1 duration and range for this model
        p_transit_model->SetTransitCellG1Duration(2.0);
        p_transit_model->SetRange(2.0);

        CustomUniformG1GenerationalCellCycleModel* p_diff_model = new CustomUniformG1GenerationalCellCycleModel;

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->InitialiseCellCycleModel();

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(
            4.0 * (p_stem_model->GetStemCellG1Duration() + p_stem_model->GetSG2MDuration()), 2 * num_steps);

        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three
            // random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 3.19525);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 2.18569);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132); // any old number
        }

        CustomUniformG1GenerationalCellCycleModel* p_hepa_one_model = new CustomUniformG1GenerationalCellCycleModel;

        // Change G1 duration and range for this model
        p_hepa_one_model->SetStemCellG1Duration(3.0);
        p_hepa_one_model->SetRange(4.0);

        CellPtr p_hepa_one_cell(new Cell(p_healthy_state, p_hepa_one_model));
        p_hepa_one_cell->SetCellProliferativeType(p_stem_type);
        p_hepa_one_cell->InitialiseCellCycleModel();


        for (unsigned i = 0; i < num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            CheckReadyToDivideAndPhaseIsUpdated(p_hepa_one_model, 3.86076);
        }
    }

    void TestArchiveCustomUniformG1GenerationalCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CustomUniformG1GenerationalCellCycleModel.arch";

        // We will also test that the random number generator is archived correctly
        double random_number_test = 0.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new CustomUniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);
            static_cast<CustomUniformG1GenerationalCellCycleModel*>(p_model)->SetTransitCellG1Duration(1.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            delete p_model;
            SimulationTime::Destroy();

            random_number_test = RandomNumberGenerator::Instance()->ranf();
            RandomNumberGenerator::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-6);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestCellCycleModelOutputParameters()
    {
        // TODO Look at TestPerimeterRatioWriter
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with CustomUniformG1GenerationalCellCycleModel
        CustomUniformG1GenerationalCellCycleModel custom_uniform_distributed_generation_based_cell_cycle_model;
        TS_ASSERT_EQUALS(custom_uniform_distributed_generation_based_cell_cycle_model.GetIdentifier(), "CustomUniformG1GenerationalCellCycleModel");

        out_stream custom_uniform_distributed_generation_based_parameter_file = output_file_handler.OpenOutputFile("custom_uniform_distributed_generation_based_results.parameters");
        custom_uniform_distributed_generation_based_cell_cycle_model.OutputCellCycleModelParameters(custom_uniform_distributed_generation_based_parameter_file);
        custom_uniform_distributed_generation_based_parameter_file->close();

        // Generate expected test file
        out_stream test_file =
            output_file_handler.OpenOutputFile("test_custom_uniform_distributed_generation_based_results.parameters");

        *test_file << "\t\t\t<Range>" << 2.0 << "</Range>\n";
        *test_file << "\t\t\t<MaxTransitGenerations>" << 3 << "</MaxTransitGenerations>\n";
        *test_file << "\t\t\t<StemCellG1Duration>" << 14.0 << "</StemCellG1Duration>\n";
        *test_file << "\t\t\t<TransitCellG1Duration>" << 2.0 << "</TransitCellG1Duration>\n";
        *test_file << "\t\t\t<SDuration>" << 5.0 << "</SDuration>\n";
        *test_file << "\t\t\t<G2Duration>" << 4.0 << "</G2Duration>\n";
        *test_file << "\t\t\t<MDuration>" << 1.0 << "</MDuration>\n";

        test_file->close();


        {
            FileFinder generated_file = output_file_handler.FindFile("custom_uniform_distributed_generation_based_results.parameters");
            FileFinder reference_file = output_file_handler.FindFile("test_custom_uniform_distributed_generation_based_results.parameters");
            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};
#endif // TESTCUSTOMUNIFORMG1GENERATIONALCELLCYCLEMODEL_HPP
