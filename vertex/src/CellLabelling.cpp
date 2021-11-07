#include "CellsGenerator.hpp"
#include "CellAncestor.hpp"

#include "CellLabelling.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"

void RandomlyLabelCells(std::list<CellPtr>& rCells,
        boost::shared_ptr<AbstractCellProperty> pLabel, double BCellsRatio,
        double b_cell_target_area, double b_cell_g1_duration,
        double b_cell_g2_duration, unsigned b_cell_max_generations)
{
    // Check that rCells is not empty (RandomNumberGenerator::Shuffle doesn't
    // work with zero size).
    assert(rCells.size() > 0);

    // Generate random permutation
    std::vector<unsigned> shuffled_cells;
    RandomNumberGenerator::Instance()->Shuffle(rCells.size(), shuffled_cells);

    unsigned cutoff = BCellsRatio * rCells.size();
    unsigned idx = 0;

    for (auto cell_iter = rCells.begin();
         cell_iter != rCells.end();
         ++cell_iter)
    {
        if (shuffled_cells[idx] < cutoff)
        {
            // Set ancestor
            MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (1)); // cell type B is ancestor 1
            (*cell_iter)->SetAncestor(p_cell_ancestor);

           (*cell_iter)->AddCellProperty(pLabel);
           (*cell_iter)->GetCellData()->SetItem("target area", b_cell_target_area);

            ExponentialG1GenerationalCellCycleModel* p_cell_cycle_model =
                static_cast<ExponentialG1GenerationalCellCycleModel*>(
                        (*cell_iter)->GetCellCycleModel());

            // Set G1 duration (exponentially distributed)
            p_cell_cycle_model->SetTransitCellG1Duration(b_cell_g1_duration);

            // set M and S duration effectively to 0
            p_cell_cycle_model->SetMDuration(1e-12);
            p_cell_cycle_model->SetSDuration(1e-12);

            // Set G2 duration
            p_cell_cycle_model->SetG2Duration(b_cell_g2_duration);

            // Set maximum number of generations
            p_cell_cycle_model->SetMaxTransitGenerations(
                    b_cell_max_generations);

            // Set birth time
            double birth_time = -
                p_cell_cycle_model->GetAverageTransitCellCycleTime() *
                RandomNumberGenerator::Instance()->ranf();
            (*cell_iter)->SetBirthTime(birth_time);

            p_cell_cycle_model->Initialise();
        }

        idx++;
    }
}

void DeterministicallyLabelCells(std::list<CellPtr>& rCells,
        boost::shared_ptr<AbstractCellProperty> pLabel, double BCellsRatio,
        double b_cell_target_area, double b_cell_g1_duration,
        double b_cell_g2_duration, unsigned b_cell_max_generations)
{
    unsigned cutoff = BCellsRatio * rCells.size();
    unsigned idx = 0;

    for (auto cell_iter = rCells.begin();
         cell_iter != rCells.end();
         ++cell_iter)
    {
        if (idx == cutoff)
            break;

        // Set ancestor
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (1)); // cell type B is ancestor 1
        (*cell_iter)->SetAncestor(p_cell_ancestor);

        (*cell_iter)->AddCellProperty(pLabel);
        (*cell_iter)->GetCellData()->SetItem("target area", b_cell_target_area);

        ExponentialG1GenerationalCellCycleModel* p_cell_cycle_model =
            static_cast<ExponentialG1GenerationalCellCycleModel*>(
                    (*cell_iter)->GetCellCycleModel());

        // Set G1 duration (exponentially distributed)
        p_cell_cycle_model->SetTransitCellG1Duration(b_cell_g1_duration);

        // set M and S duration effectively to 0
        p_cell_cycle_model->SetMDuration(1e-12);
        p_cell_cycle_model->SetSDuration(1e-12);

        // Set G2 duration
        p_cell_cycle_model->SetG2Duration(b_cell_g2_duration);

        // Set maximum number of generations
        p_cell_cycle_model->SetMaxTransitGenerations(
                b_cell_max_generations);

        // Set birth time
        double birth_time = -
            p_cell_cycle_model->GetAverageTransitCellCycleTime() *
            RandomNumberGenerator::Instance()->ranf();
        (*cell_iter)->SetBirthTime(birth_time);

        p_cell_cycle_model->Initialise();

        idx++;
    }
}
