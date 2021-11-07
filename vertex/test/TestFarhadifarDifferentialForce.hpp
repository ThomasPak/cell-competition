#ifndef TESTFARHADIFARDIFFERENTIALFORCE_HPP_
#define TESTFARHADIFARDIFFERENTIALFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"

#include "SmartPointers.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"

#include "FarhadifarDifferentialForce.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestFarhadifarDifferentialForce : public AbstractCellBasedTestSuite
{
public:

    void TestFarhadifarDifferentialForceMethods1()
    {
        // This is the same test as for other vertex based forces. It comprises
        // a sanity check that forces point in the right direction.
        // Construct a 2D vertex mesh consisting of a single element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        std::vector<double> angles = std::vector<double>(num_nodes);

        for (unsigned i=0; i<num_nodes; i++)
        {
            angles[i] = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, true, cos(angles[i]), sin(angles[i])));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up the cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(-1.0);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Create a force system
        FarhadifarDifferentialForce<2> force;

        // Test get/set methods
        TS_ASSERT_DELTA(force.GetCellAreaElasticityParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellAreaElasticityParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetCellPerimeterContractilityParameter(), 0.04, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellPerimeterContractilityParameter(), 0.04, 1e-12);
        TS_ASSERT_DELTA(force.GetCellCellLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellLabelledCellLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellCellLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetCellBoundaryLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellBoundaryLineTensionParameter(), 0.12, 1e-12);

        force.SetCellAreaElasticityParameter(5.8);
        force.SetLabelledCellAreaElasticityParameter(4.6);
        force.SetCellPerimeterContractilityParameter(17.9);
        force.SetLabelledCellPerimeterContractilityParameter(16.6);
        force.SetCellCellLineTensionParameter(0.5);
        force.SetLabelledCellLabelledCellLineTensionParameter(0.7);
        force.SetLabelledCellCellLineTensionParameter(0.4);
        force.SetCellBoundaryLineTensionParameter(0.6);
        force.SetLabelledCellBoundaryLineTensionParameter(0.3);

        TS_ASSERT_DELTA(force.GetCellAreaElasticityParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellAreaElasticityParameter(), 4.6, 1e-12);
        TS_ASSERT_DELTA(force.GetCellPerimeterContractilityParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellPerimeterContractilityParameter(), 16.6, 1e-12);
        TS_ASSERT_DELTA(force.GetCellCellLineTensionParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellLabelledCellLineTensionParameter(), 0.7, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellCellLineTensionParameter(), 0.4, 1e-12);
        TS_ASSERT_DELTA(force.GetCellBoundaryLineTensionParameter(), 0.6, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellBoundaryLineTensionParameter(), 0.3, 1e-12);

        force.SetCellAreaElasticityParameter(1.0);
        force.SetLabelledCellAreaElasticityParameter(1.0);
        force.SetCellPerimeterContractilityParameter(0.04);
        force.SetLabelledCellPerimeterContractilityParameter(0.04);
        force.SetCellCellLineTensionParameter(0.12);
        force.SetLabelledCellLabelledCellLineTensionParameter(0.12);
        force.SetLabelledCellCellLineTensionParameter(0.12);
        force.SetCellBoundaryLineTensionParameter(0.12);
        force.SetLabelledCellBoundaryLineTensionParameter(0.12);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Currently, the Farhadifar force only works if used together with a
        // target area growth modifier
        // This tests that a meaningful error appears if we don't use a growth
        // modifier
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(cell_population),
                "You need to add an AbstractTargetAreaModifier to the simulation in order to use a FarhadifarDifferentialForce");

        // create our modifier, which sets the target areas for the cell population

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should be radially inward, with the same
        // magnitude for all nodes
        double force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.25, 2);

        // Set the cell to be necrotic
        cell_population.GetCellUsingLocationIndex(0)->StartApoptosis();

        // Reset force vector
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Increment time
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_DELTA(cell_population.GetCellUsingLocationIndex(0)->GetTimeUntilDeath(), 0.125, 1e-6);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // Now the forces should be affected
        double apoptotic_force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());
        TS_ASSERT_LESS_THAN(force_magnitude, apoptotic_force_magnitude);
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), apoptotic_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -apoptotic_force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -apoptotic_force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestFarhadifarDifferentialForceMethods2()
    {
        // This is the same test as TestFarhadifarDifferentialForceMethods1,
        // but with a labelled cell.
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 20;
        std::vector<double> angles = std::vector<double>(num_nodes);

        for (unsigned i=0; i<num_nodes; i++)
        {
            angles[i] = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, true, cos(angles[i]), sin(angles[i])));
        }

        std::vector<VertexElement<2,2>*> elements;
        elements.push_back(new VertexElement<2,2>(0, nodes));

        double cell_swap_threshold = 0.01;
        double edge_division_threshold = 2.0;
        MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, edge_division_threshold);

        // Set up the cell
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        boost::shared_ptr<AbstractCellProperty>
            p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->SetBirthTime(-1.0);
        p_cell->AddCellProperty(p_label);
        cells.push_back(p_cell);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.InitialiseCells();

        // Create a force system
        FarhadifarDifferentialForce<2> force;

        // Test get/set methods
        TS_ASSERT_DELTA(force.GetCellAreaElasticityParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellAreaElasticityParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(force.GetCellPerimeterContractilityParameter(), 0.04, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellPerimeterContractilityParameter(), 0.04, 1e-12);
        TS_ASSERT_DELTA(force.GetCellCellLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellLabelledCellLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellCellLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetCellBoundaryLineTensionParameter(), 0.12, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellBoundaryLineTensionParameter(), 0.12, 1e-12);

        force.SetCellAreaElasticityParameter(5.8);
        force.SetLabelledCellAreaElasticityParameter(4.6);
        force.SetCellPerimeterContractilityParameter(17.9);
        force.SetLabelledCellPerimeterContractilityParameter(16.6);
        force.SetCellCellLineTensionParameter(0.5);
        force.SetLabelledCellLabelledCellLineTensionParameter(0.7);
        force.SetLabelledCellCellLineTensionParameter(0.4);
        force.SetCellBoundaryLineTensionParameter(0.6);
        force.SetLabelledCellBoundaryLineTensionParameter(0.3);

        TS_ASSERT_DELTA(force.GetCellAreaElasticityParameter(), 5.8, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellAreaElasticityParameter(), 4.6, 1e-12);
        TS_ASSERT_DELTA(force.GetCellPerimeterContractilityParameter(), 17.9, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellPerimeterContractilityParameter(), 16.6, 1e-12);
        TS_ASSERT_DELTA(force.GetCellCellLineTensionParameter(), 0.5, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellLabelledCellLineTensionParameter(), 0.7, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellCellLineTensionParameter(), 0.4, 1e-12);
        TS_ASSERT_DELTA(force.GetCellBoundaryLineTensionParameter(), 0.6, 1e-12);
        TS_ASSERT_DELTA(force.GetLabelledCellBoundaryLineTensionParameter(), 0.3, 1e-12);

        force.SetCellAreaElasticityParameter(1.0);
        force.SetLabelledCellAreaElasticityParameter(1.0);
        force.SetCellPerimeterContractilityParameter(0.04);
        force.SetLabelledCellPerimeterContractilityParameter(0.04);
        force.SetCellCellLineTensionParameter(0.12);
        force.SetLabelledCellLabelledCellLineTensionParameter(0.12);
        force.SetLabelledCellCellLineTensionParameter(0.12);
        force.SetCellBoundaryLineTensionParameter(0.12);
        force.SetLabelledCellBoundaryLineTensionParameter(0.12);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Currently, the Farhadifar force only works if used together with a
        // target area growth modifier
        // This tests that a meaningful error appears if we don't use a growth
        // modifier
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(cell_population),
                "You need to add an AbstractTargetAreaModifier to the simulation in order to use a FarhadifarDifferentialForce");

        // create our modifier, which sets the target areas for the cell population

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should be radially inward, with the same
        // magnitude for all nodes
        double force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());

        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Set up simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(0.25, 2);

        // Set the cell to be necrotic
        cell_population.GetCellUsingLocationIndex(0)->StartApoptosis();

        // Reset force vector
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // The force on each node should not yet be affected by setting the cell to be apoptotic
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -force_magnitude*sin(angles[i]), 1e-4);
        }

        // Increment time
        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_DELTA(cell_population.GetCellUsingLocationIndex(0)->GetTimeUntilDeath(), 0.125, 1e-6);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }
        //#2488 workaround
        p_growth_modifier->UpdateTargetAreas(cell_population);

        force.AddForceContribution(cell_population);

        // Now the forces should be affected
        double apoptotic_force_magnitude = norm_2(cell_population.GetNode(0)->rGetAppliedForce());
        TS_ASSERT_LESS_THAN(force_magnitude, apoptotic_force_magnitude);
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_DELTA(norm_2(cell_population.GetNode(i)->rGetAppliedForce()), apoptotic_force_magnitude, 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], -apoptotic_force_magnitude*cos(angles[i]), 1e-4);
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], -apoptotic_force_magnitude*sin(angles[i]), 1e-4);
        }
    }

    void TestFarhadifarDifferentialForceTerms1()
    {
        /**
         * Here we test that the forces are applied correctly to individual nodes.
         * We apply the force to something like this:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns
        // mature cell target areas
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
        }

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        // Now let's make a FarhadifarDifferentialForce and apply it to the population
        FarhadifarDifferentialForce<2> force;

        force.AddForceContribution(cell_population);

        c_vector<double, 2> applied_force_0;
        applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1;
        applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // If this is a Farhadifar force, this will be the force at the vertices
        TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[1], 6.76, 1e-10);
    }

    void TestFarhadifarDifferentialForceTerms2()
    {
        // This is the same test as TestFarhadifarDifferentialForceTerms1,
        // but with two labelled cells.
        /**
         * Here we test that the forces are applied correctly to individual nodes.
         * We apply the force to something like this:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns
        // mature cell target areas
        // Also label cells
        boost::shared_ptr<AbstractCellProperty>
            p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
            cell_iter->AddCellProperty(p_label);
            cell_iter->GetCellData()->SetItem("target area", 1.5);
        }

        // Now let's make a FarhadifarDifferentialForce and apply it to the population
        FarhadifarDifferentialForce<2> force;
        force.SetLabelledCellAreaElasticityParameter(1.5);
        force.SetLabelledCellLabelledCellLineTensionParameter(0.06);
        force.SetLabelledCellCellLineTensionParameter(0.14);
        force.SetLabelledCellBoundaryLineTensionParameter(0.10);
        force.SetLabelledCellPerimeterContractilityParameter(0.11);

        force.AddForceContribution(cell_population);

        c_vector<double, 2> applied_force_0;
        applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1;
        applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // If this is a Farhadifar force, this will be the force at the vertices
        TS_ASSERT_DELTA(applied_force_0[0], 4.73, 1e-10);
        TS_ASSERT_DELTA(applied_force_0[1], 4.73, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[1], 9.32, 1e-10);
    }

    void TestFarhadifarDifferentialForceTerms3()
    {
        // This is the same test as TestFarhadifarDifferentialForceTerms1,
        // but with one unlabelled cells and one labelled
        /**
         * Here we test that the forces are applied correctly to individual nodes.
         * We apply the force to something like this:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns
        // mature cell target areas
        // Also label cells
        boost::shared_ptr<AbstractCellProperty>
            p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
            cell_iter->GetCellData()->SetItem("target area", 1.0);
        }

        // Set labelled cell properties
        auto cell_iter = cell_population.Begin();
        ++cell_iter;
        cell_iter->AddCellProperty(p_label);
        cell_iter->GetCellData()->SetItem("target area", 1.5);

        // Now let's make a FarhadifarDifferentialForce and apply it to the population
        FarhadifarDifferentialForce<2> force;
        force.SetLabelledCellAreaElasticityParameter(1.5);
        force.SetLabelledCellLabelledCellLineTensionParameter(0.06);
        force.SetLabelledCellCellLineTensionParameter(0.14);
        force.SetLabelledCellBoundaryLineTensionParameter(0.10);
        force.SetLabelledCellPerimeterContractilityParameter(0.11);

        force.AddForceContribution(cell_population);

        c_vector<double, 2> applied_force_0;
        applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1;
        applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // If this is a Farhadifar force, this will be the force at the vertices
        TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[0], 1.29, 1e-10);
        TS_ASSERT_DELTA(applied_force_1[1], 8.09, 1e-10);
    }

    void TestFarhadifarDifferentialForceInSimulation1()
    {
        /**
         * This is the same test as TestFarhadifarDifferentialForceTerms1, just
         * that now we don't check that the applied forces are calculated
         * correctly, but rather that in a simulation the displacement of
         * vertices is as we expect.
         *
         * This is the mesh:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
        }

        MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
        p_growth_modifier->UpdateTargetAreas(cell_population);

        // Now let's make a FarhadifarDifferentialForce and add it to the simulation.
        MAKE_PTR(FarhadifarDifferentialForce<2>, p_force);

        // We need to reset the cell rearrangement threshold - vertex movements are kept below that threshold
        cell_population.rGetMesh().SetCellRearrangementThreshold(0.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFarhadifarDifferentialForce");
        simulator.SetEndTime(0.01);
        simulator.SetDt(0.01);
        simulator.AddForce(p_force);

        simulator.Solve();

        c_vector<double, 2> applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // New Location = Old Location + (Dt * applied force), since viscosity should be one
        c_vector<double, 2> expected_new_node_location_0;
        expected_new_node_location_0[0] = 0.0+0.01*3.44;
        expected_new_node_location_0[1] = 0.0+0.01*3.44;
        c_vector<double, 2> expected_new_node_location_1;
        expected_new_node_location_1[0] = 2.0 + 0.01*0.0;
        expected_new_node_location_1[1] = 0.0 + 0.01*6.76;

        // If this is a Farhadifar force, this will be the location of the first two vertices.
        TS_ASSERT_DELTA(expected_new_node_location_0[0], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_0[1], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[1], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[0], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[1], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[1], 1e-10);

    }

    void TestFarhadifarDifferentialForceInSimulation2()
    {
        /**
         * This is the same test as
         * TestFarhadifarDifferentialForceInSimulation1, but with two labelled
         * cells.
         *
         * This is the mesh:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns
        // mature cell target areas
        boost::shared_ptr<AbstractCellProperty>
            p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
            cell_iter->AddCellProperty(p_label);
            cell_iter->GetCellData()->SetItem("target area", 1.5);
        }

        // Now let's make a FarhadifarDifferentialForce and add it to the simulation.
        MAKE_PTR(FarhadifarDifferentialForce<2>, p_force);
        p_force->SetLabelledCellAreaElasticityParameter(1.5);
        p_force->SetLabelledCellLabelledCellLineTensionParameter(0.06);
        p_force->SetLabelledCellCellLineTensionParameter(0.14);
        p_force->SetLabelledCellBoundaryLineTensionParameter(0.10);
        p_force->SetLabelledCellPerimeterContractilityParameter(0.11);

        // We need to reset the cell rearrangement threshold - vertex movements are kept below that threshold
        cell_population.rGetMesh().SetCellRearrangementThreshold(0.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFarhadifarDifferentialForce");
        simulator.SetEndTime(0.01);
        simulator.SetDt(0.01);
        simulator.AddForce(p_force);

        simulator.Solve();

        c_vector<double, 2> applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // New Location = Old Location + (Dt * applied force), since viscosity should be one
        c_vector<double, 2> expected_new_node_location_0;
        expected_new_node_location_0[0] = 0.0+0.01*4.73;
        expected_new_node_location_0[1] = 0.0+0.01*4.73;
        c_vector<double, 2> expected_new_node_location_1;
        expected_new_node_location_1[0] = 2.0 + 0.01*0.0;
        expected_new_node_location_1[1] = 0.0 + 0.01*9.32;

        // If this is a Farhadifar force, this will be the location of the first two vertices.
        TS_ASSERT_DELTA(expected_new_node_location_0[0], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_0[1], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[1], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[0], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[1], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[1], 1e-10);

    }

    void TestFarhadifarDifferentialForceInSimulation3()
    {
        /**
         * This is the same test as
         * TestFarhadifarDifferentialForceInSimulation1, but with one
         * unlabelled cell and one labelled cell.
         *
         * This is the mesh:
         *  . ____ . ____ .
         *  |      |      |
         *  |      |      |
         *  . ____ . ____ .
         */
        std::vector<Node<2>*> nodes;
        // the boolean says wether the node is a boundary node or not
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
        nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
        nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

        // make two square elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
        unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

        for (unsigned i=0; i<4; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
        }

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

        // Get a cell population
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        // Set the birth time to -5 such that the target area modifier assigns
        // mature cell target areas
        boost::shared_ptr<AbstractCellProperty>
            p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            cell_iter->SetBirthTime(-5.0);
            cell_iter->GetCellData()->SetItem("target area", 1.0);
        }

        // Set labelled cell properties
        auto cell_iter = cell_population.Begin();
        ++cell_iter;
        cell_iter->AddCellProperty(p_label);
        cell_iter->GetCellData()->SetItem("target area", 1.5);

        // Now let's make a FarhadifarDifferentialForce and add it to the simulation.
        MAKE_PTR(FarhadifarDifferentialForce<2>, p_force);
        p_force->SetLabelledCellAreaElasticityParameter(1.5);
        p_force->SetLabelledCellLabelledCellLineTensionParameter(0.06);
        p_force->SetLabelledCellCellLineTensionParameter(0.14);
        p_force->SetLabelledCellBoundaryLineTensionParameter(0.10);
        p_force->SetLabelledCellPerimeterContractilityParameter(0.11);

        // We need to reset the cell rearrangement threshold - vertex movements are kept below that threshold
        cell_population.rGetMesh().SetCellRearrangementThreshold(0.5);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFarhadifarDifferentialForce");
        simulator.SetEndTime(0.01);
        simulator.SetDt(0.01);
        simulator.AddForce(p_force);

        simulator.Solve();

        c_vector<double, 2> applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        c_vector<double, 2> applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        // New Location = Old Location + (Dt * applied force), since viscosity should be one
        c_vector<double, 2> expected_new_node_location_0;
        expected_new_node_location_0[0] = 0.0+0.01*3.44;
        expected_new_node_location_0[1] = 0.0+0.01*3.44;
        c_vector<double, 2> expected_new_node_location_1;
        expected_new_node_location_1[0] = 2.0 + 0.01*1.29;
        expected_new_node_location_1[1] = 0.0 + 0.01*8.09;

        // If this is a Farhadifar force, this will be the location of the first two vertices.
        TS_ASSERT_DELTA(expected_new_node_location_0[0], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_0[1], (cell_population.rGetMesh().GetNode(0)->rGetLocation())[1], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[0], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[0], 1e-10);
        TS_ASSERT_DELTA(expected_new_node_location_1[1], (cell_population.rGetMesh().GetNode(1)->rGetLocation())[1], 1e-10);

    }
};

#endif // TESTFARHADIFARDIFFERENTIALFORCE_HPP
