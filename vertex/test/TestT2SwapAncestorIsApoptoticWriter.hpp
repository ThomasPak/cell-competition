#ifndef TESTT2SWAPANCESTORISAPOPTOTICWRITER_HPP
#define TESTT2SWAPANCESTORISAPOPTOTICWRITER_HPP

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

// Cell population writers
#include "T2SwapAncestorIsApoptoticWriter.hpp"

// Files to create populations
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CustomVertexBasedCellPopulation.hpp"
#include "CustomT2SwapCellKiller.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestT2SwapAncestorIsApoptoticWriter : public AbstractCellBasedTestSuite
{
public:

    void TestT2SwapAncestorIsApoptoticWriter1()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        CustomVertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestT2SwapAncestorIsApoptoticWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a T2SwapAncestorIsApoptoticWriter and test that the correct output is generated
        T2SwapAncestorIsApoptoticWriter<2,2> t2_swaps_writer;
        t2_swaps_writer.OpenOutputFile(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2AncestorsIsApoptotics.dat", "cell_based/test/data/TestCellPopulationWriters/T2SwapLocations.dat").CompareFiles();

        // Test that we can append to files
        t2_swaps_writer.OpenOutputFileForAppend(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2AncestorsIsApoptotics.dat", "cell_based/test/data/TestCellPopulationWriters/T2SwapLocations_twice.dat").CompareFiles();

        {
            // Coverage of the Visit() method when called on a MeshBasedCellPopulation
            HoneycombMeshGenerator tet_generator(5, 5, 0);
            MutableMesh<2,2>* p_tet_mesh = tet_generator.GetMesh();
            std::vector<CellPtr> mesh_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
            mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
            MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&mesh_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a CaBasedCellPopulation
            PottsMeshGenerator<2> ca_based_generator(5, 0, 0, 5, 0, 0);
            PottsMesh<2>* p_ca_based_mesh = ca_based_generator.GetMesh();
            std::vector<CellPtr> ca_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> ca_based_cells_generator;
            ca_based_cells_generator.GenerateBasic(ca_based_cells, 5);
            std::vector<unsigned> location_indices;
            location_indices.push_back(7);
            location_indices.push_back(11);
            location_indices.push_back(12);
            location_indices.push_back(13);
            location_indices.push_back(17);
            CaBasedCellPopulation<2> ca_based_cell_population(*p_ca_based_mesh, ca_based_cells, location_indices);

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&ca_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a NodeBasedCellPopulation
            std::vector<Node<2>* > node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&node_based_cell_population));

            // Tidy up
            delete node_based_nodes[0];
            delete node_based_nodes[1];
        }

        {
            // Coverage of the Visit() method when called on a PottsBasedCellPopulation
            PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
            PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
            std::vector<CellPtr> potts_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
            potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
            PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&potts_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a VertexBasedCellPopulation
            HoneycombVertexMeshGenerator vertex_generator(4, 6);
            MutableVertexMesh<2,2>* p_vertex_mesh = vertex_generator.GetMesh();
            std::vector<CellPtr> vertex_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> vertex_based_cells_generator;
            vertex_based_cells_generator.GenerateBasic(vertex_based_cells,
                    p_vertex_mesh->GetNumElements());
            VertexBasedCellPopulation<2>
                vertex_based_cell_population(*p_vertex_mesh, vertex_based_cells);

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&vertex_based_cell_population));
        }
    }

    void TestT2SwapAncestorIsApoptoticWriter2()
    {
        EXIT_IF_PARALLEL;

        // Make 6 nodes to assign to four elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 0.5, 0.5));
        nodes.push_back(new Node<2>(3, false, 0.1, 0.05));
        nodes.push_back(new Node<2>(4, false, 0.9, 0.05));
        nodes.push_back(new Node<2>(5, false, 0.5, 0.475));

        // Make one triangular and three trapezium elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3;
        unsigned node_indices_elem_0[3] = {3, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 5, 4};
        unsigned node_indices_elem_2[4] = {2, 0, 3, 5};
        unsigned node_indices_elem_3[4] = {0, 1, 4, 3};
        for (unsigned i=0; i<4; i++)
        {
            if (i < 3)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            }
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
            nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        vertex_mesh.SetT2Threshold(0.01);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        CustomVertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set cell ancestor to location indices
        cell_population.SetCellAncestorsToLocationIndices();

        // Call apoptosis on T2 cell
        auto cell_iter = cell_population.Begin();
        cell_iter->StartApoptosis(false);

        // The population should have 4 cells
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u)

        // Give the Population to the cell killer
        CustomT2SwapCellKiller<2> cell_killer(&cell_population);

        // We move the inner vertices a bit inwards
        c_vector<double, 2>& new_location_0 = vertex_elements[0]->GetNode(0)->rGetModifiableLocation();
        new_location_0(0) = 0.499;
        new_location_0(1) = 0.249;

        c_vector<double, 2>& new_location_1 = vertex_elements[0]->GetNode(1)->rGetModifiableLocation();
        new_location_1(0) = 0.501;
        new_location_1(1) = 0.249;

        c_vector<double, 2>& new_location_2 = vertex_elements[0]->GetNode(2)->rGetModifiableLocation();
        new_location_2(0) = 0.5;
        new_location_2(1) = 0.251;

        // Perform swaps
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Create an output directory for the writer
        std::string output_directory = "TestT2SwapAncestorIsApoptoticWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create test files
        out_stream test_file_1 =
            output_file_handler.OpenOutputFile("test_ancestors_1.dat");
        *test_file_1 << "0\t1\t0\tapoptosis\t\n";
        test_file_1->close();

        out_stream test_file_2 =
            output_file_handler.OpenOutputFile("test_ancestors_2.dat");
        *test_file_2 << "0\t1\t0\tapoptosis\t\n0\t0\t\n";
        test_file_2->close();

        // Create a T2SwapAncestorIsApoptoticWriter and test that the correct output is generated
        T2SwapAncestorIsApoptoticWriter<2,2> t2_swaps_writer;
        t2_swaps_writer.OpenOutputFile(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2AncestorsIsApoptotics.dat", results_dir
                + "test_ancestors_1.dat").CompareFiles();

        // Test that we can append to files
        t2_swaps_writer.OpenOutputFileForAppend(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2AncestorsIsApoptotics.dat", results_dir
                + "test_ancestors_2.dat").CompareFiles();
    }

    void TestT2SwapAncestorIsApoptoticWriter3()
    {
        EXIT_IF_PARALLEL;

        /**
         * Create a mesh comprising ten nodes contained in six elements, two of which are small triangles,
         * as shown below. We will test that the CheckAndLabelCellsForApoptosisOrDeath() method works correctly in the case where multiple
         * T2 swaps are required. After remeshing, the two triangular elements should be removed from the mesh.
         *
         *   ________          _________
         *  |\      /|        |\       /|
         *  | |\__/| | -----> | \_____/ |
         *  | |/  \| |        | /     \ |
         *  |/______\|        |/_______\|
         *
         */
        std::vector<Node<2>*> nodes;
        // Create the nodes
        // The boolean is true for boundary nodes
        nodes.push_back(new Node<2>(0, true,  0.0, 0.0));
        nodes.push_back(new Node<2>(1, true,  1.0, 0.0));
        nodes.push_back(new Node<2>(2, true,  1.0, 1.0));
        nodes.push_back(new Node<2>(3, true,  0.0, 1.0));
        nodes.push_back(new Node<2>(4, false, 0.3, 0.45));
        nodes.push_back(new Node<2>(5, false, 0.4, 0.5));
        nodes.push_back(new Node<2>(6, false, 0.3, 0.55));
        nodes.push_back(new Node<2>(7, false, 0.6, 0.5));
        nodes.push_back(new Node<2>(8, false, 0.7, 0.45));
        nodes.push_back(new Node<2>(9, false, 0.7, 0.55));

        // Construct elements from the nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1, nodes_elem_2, nodes_elem_3, nodes_elem_4, nodes_elem_5;
        unsigned node_indices_elem_0[6] = {0, 1, 8, 7, 5, 4};
        unsigned node_indices_elem_1[6] = {2, 3, 6, 5, 7, 9};
        unsigned node_indices_elem_2[4] = {1, 2, 9, 8};
        unsigned node_indices_elem_3[4] = {0, 4, 6, 3};
        unsigned node_indices_elem_4[3] = {4, 5, 6};
        unsigned node_indices_elem_5[3] = {7, 8, 9};

        for (unsigned i=0; i<6; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            if (i < 4)
            {
                nodes_elem_2.push_back(nodes[node_indices_elem_2[i]]);
                nodes_elem_3.push_back(nodes[node_indices_elem_3[i]]);
            }
            if (i < 3)
            {
                nodes_elem_4.push_back(nodes[node_indices_elem_4[i]]);
                nodes_elem_5.push_back(nodes[node_indices_elem_5[i]]);
            }
        }

        // create a mesh
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));
        vertex_elements.push_back(new VertexElement<2,2>(2, nodes_elem_2));
        vertex_elements.push_back(new VertexElement<2,2>(3, nodes_elem_3));
        vertex_elements.push_back(new VertexElement<2,2>(4, nodes_elem_4));
        vertex_elements.push_back(new VertexElement<2,2>(5, nodes_elem_5));

        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        // Set the remeshing threshold parameter values so that T2 swaps do occur, but T1 swaps do not
        vertex_mesh.SetT2Threshold(0.01);
        vertex_mesh.SetCellRearrangementThreshold(0.00001);

        // Test that the numbers of nodes and elements are correct
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        CustomVertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set cell ancestor to location indices
        cell_population.SetCellAncestorsToLocationIndices();

        // Call apoptosis on first cell to do T2 swap (index 4)
        auto cell_iter = cell_population.Begin();
        for (int i = 0; i < 4; i++)
        {
            ++cell_iter;
        }
        cell_iter->StartApoptosis();

        // the population should have 6 cells
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 6u)

        // Give the Population to the cell killer
        CustomT2SwapCellKiller<2> cell_killer(&cell_population);
        // Perform swaps
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Create an output directory for the writer
        std::string output_directory = "TestT2SwapAncestorIsApoptoticWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create test files
        out_stream test_file_3 =
            output_file_handler.OpenOutputFile("test_ancestors_3.dat");
        *test_file_3 << "0\t2\t4\tapoptosis\t5\textrusion\t\n";
        test_file_3->close();

        out_stream test_file_4 =
            output_file_handler.OpenOutputFile("test_ancestors_4.dat");
        *test_file_4 << "0\t2\t4\tapoptosis\t5\textrusion\t\n0\t0\t\n";
        test_file_4->close();

        // Create a T2SwapAncestorIsApoptoticWriter and test that the correct output is generated
        T2SwapAncestorIsApoptoticWriter<2,2> t2_swaps_writer;
        t2_swaps_writer.OpenOutputFile(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2AncestorsIsApoptotics.dat", results_dir
                + "test_ancestors_3.dat").CompareFiles();

        // Test that we can append to files
        t2_swaps_writer.OpenOutputFileForAppend(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2AncestorsIsApoptotics.dat", results_dir
                + "test_ancestors_4.dat").CompareFiles();
    }

    void TestT2SwapAncestorIsApoptoticWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "T2SwapAncestorIsApoptoticWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new T2SwapAncestorIsApoptoticWriter<2,2>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cell_writer;
            delete p_cell_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_writer_2;
            delete p_cell_writer_2;
       }
    }
};


#endif // TESTT2SWAPANCESTORISAPOPTOTICWRITER_HPP
