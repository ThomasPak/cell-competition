#ifndef TESTPERIMETERRATIOWRITER_HPP
#define TESTPERIMETERRATIOWRITER_HPP

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "PerimeterRatioWriter.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "UniformCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellAncestorWriter.hpp"
#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"
#include "CellAncestor.hpp"

#include "PetscSetupAndFinalize.hpp"

static constexpr double EDGE_LENGTH = pow(3.0, -0.5);

class TestPerimeterRatioWriter : public AbstractCellBasedTestSuite
{
    public:

        // try
        void TestPerimeterRatio1()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Add cell writer
            cell_population.AddCellWriter<CellAncestorWriter>();

            // Add population writer
            cell_population.AddPopulationWriter<PerimeterRatioWriter>();

            // Set cell ancestors
            cell_population.SetCellAncestorsToLocationIndices();

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            std::string output_directory = "TestPerimeterRatio1";
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(0.01);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Generate expected test file
            OutputFileHandler output_file_handler(output_directory, false);
            out_stream test_file =
                output_file_handler.OpenOutputFile("test_perimeterratio.dat");

            double t = 0.0;
            for (int i = 0; i < 6; i++)
            {
                *test_file << t << "\t";

                for (int ancestor = 0; ancestor < 9; ancestor++)
                {
                    *test_file << ancestor << " " << 0.0 << "\t";
                }

                *test_file << "\n";

                // default timestep defined in VertexBasedCellPopulation.cpp::816
                t += 0.002;
            }

            test_file->close();

            // Compare files
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
            FileComparison(results_dir + "results_from_time_0/perimeterratio.dat",
                    results_dir + "test_perimeterratio.dat").CompareFiles();
        }

        void TestPerimeterRatio2()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Add cell writer
            cell_population.AddCellWriter<CellAncestorWriter>();

            // Add population writer
            cell_population.AddPopulationWriter<PerimeterRatioWriter>();

            // Set all cells to same cell ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (0));
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->SetAncestor(p_cell_ancestor);
            }

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            std::string output_directory = "TestPerimeterRatio2";
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(0.01);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Generate expected test file
            OutputFileHandler output_file_handler(output_directory, false);
            out_stream test_file =
                output_file_handler.OpenOutputFile("test_perimeterratio.dat");

            double t = 0.0;
            double perimeter_ratio = (2.0 + 4.0 + 3.0 + 5.0 + 6.0 + 3.0 + 2.0 +
                    4.0 + 3.0) / 6.0 / 9.0;
            for (int i = 0; i < 6; i++)
            {
                *test_file << t << "\t";
                *test_file << 0 << " " << perimeter_ratio << "\t";
                *test_file << "\n";

                t += 0.002;
            }

            test_file->close();

            // Compare files
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
            FileComparison(results_dir + "results_from_time_0/perimeterratio.dat",
                    results_dir + "test_perimeterratio.dat").CompareFiles();
        }

        void TestPerimeterRatio3()
        {
            EXIT_IF_PARALLEL;

            // Create a regular vertex mesh
            HoneycombVertexMeshGenerator generator(3, 3);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Add cell writer
            cell_population.AddCellWriter<CellAncestorWriter>();

            // Add population writer
            cell_population.AddPopulationWriter<PerimeterRatioWriter>();

            // Set all cells to same cell ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor_a, (0));
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor_b, (1));
            unsigned counter = 0;
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                if (counter != 4)
                {
                    cell_iter->SetAncestor(p_cell_ancestor_a);
                }
                else
                {
                    cell_iter->SetAncestor(p_cell_ancestor_b);
                }
                ++counter;
            }

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            std::string output_directory = "TestPerimeterRatio3";
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(0.01);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Generate expected test file
            OutputFileHandler output_file_handler(output_directory, false);
            out_stream test_file =
                output_file_handler.OpenOutputFile("test_perimeterratio.dat");

            double t = 0.0;
            double perimeter_ratio_a = (2.0 + 4.0 + 3.0 + 5.0 + 6.0 + 3.0 + 2.0 +
                    4.0 + 3.0 - 12.0) / 6.0 / 8.0;
            double perimeter_ratio_b = 0.0;
            for (int i = 0; i < 6; i++)
            {
                *test_file << t << "\t";
                *test_file << 0 << " " << perimeter_ratio_a << "\t";
                *test_file << 1 << " " << perimeter_ratio_b << "\t";
                *test_file << "\n";

                t += 0.002;
            }

            test_file->close();

            // Compare files
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
            FileComparison(results_dir + "results_from_time_0/perimeterratio.dat",
                    results_dir + "test_perimeterratio.dat").CompareFiles();
        }

        void TestPerimeterRatio4()
        {
            EXIT_IF_PARALLEL;

            // Create custom mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
            nodes.push_back(new Node<2>(2, true, 3.0, 0.0));
            nodes.push_back(new Node<2>(3, true, 0.0, 1.0));
            nodes.push_back(new Node<2>(4, true, 2.0, 1.0));
            nodes.push_back(new Node<2>(5, true, 3.0, 1.0));

            std::vector<VertexElement<2,2>*> elements;
            elements.push_back(new VertexElement<2,2>(0, { nodes[0], nodes[1],
                        nodes[4], nodes[3] }));
            elements.push_back(new VertexElement<2,2>(1, { nodes[1], nodes[2],
                        nodes[5], nodes[4] }));

            double cell_swap_threshold = 0.01;
            double edge_division_threshold = 2.0;
            MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold,
                    edge_division_threshold);

            // Create some cells with a differentiated cell type so that they do not divide.
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_state);
            cells_generator.GenerateBasic(cells, mesh.GetNumElements(),
                    std::vector<unsigned>(), p_state);

            // Create cell-based population object
            VertexBasedCellPopulation<2> cell_population(mesh, cells);

            // Add cell writer
            cell_population.AddCellWriter<CellAncestorWriter>();

            // Add population writer
            cell_population.AddPopulationWriter<PerimeterRatioWriter>();

            // Set all cells to same cell ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (0));
            for (auto cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->SetAncestor(p_cell_ancestor);
            }

            // Create and configure cell-based simulation
            OffLatticeSimulation<2> simulator(cell_population);
            std::string output_directory = "TestPerimeterRatio4";
            simulator.SetOutputDirectory(output_directory);
            simulator.SetEndTime(0.01);

            // Run simulation
            TS_ASSERT_THROWS_NOTHING(simulator.Solve());

            // Generate expected test file
            OutputFileHandler output_file_handler(output_directory, false);
            out_stream test_file =
                output_file_handler.OpenOutputFile("test_perimeterratio.dat");

            double t = 0.0;
            double perimeter_ratio = 7.0 / 36.0;
            for (int i = 0; i < 6; i++)
            {
                *test_file << t << "\t";
                *test_file << 0 << " " << perimeter_ratio << "\t";
                *test_file << "\n";

                t += 0.002;
            }

            test_file->close();

            // Compare files
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
            FileComparison(results_dir + "results_from_time_0/perimeterratio.dat",
                    results_dir + "test_perimeterratio.dat").CompareFiles();
        }
};

#endif // TESTPERIMETERRATIOWRITER_HPP
