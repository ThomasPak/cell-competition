#ifndef CELLLABELLING_HPP_
#define CELLLABELLING_HPP_

#include "CellLabel.hpp"
#include "SmartPointers.hpp"

void RandomlyLabelCells(std::list<CellPtr>& rCells,
        boost::shared_ptr<AbstractCellProperty> pLabel, double BCellsRatio,
        double b_cell_target_area, double b_cell_g1_duration,
        double b_cell_g2_duration, unsigned b_cell_max_generations);

void DeterministicallyLabelCells(std::list<CellPtr>& rCells,
        boost::shared_ptr<AbstractCellProperty> pLabel, double BCellsRatio,
        double b_cell_target_area, double b_cell_g1_duration,
        double b_cell_g2_duration, unsigned b_cell_max_generations);

#endif // CELLLABELLING_HPP_
