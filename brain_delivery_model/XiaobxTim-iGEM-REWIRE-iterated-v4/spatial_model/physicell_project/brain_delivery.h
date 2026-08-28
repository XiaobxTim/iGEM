#pragma once

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

void setup_microenvironment();
void create_cell_types();
void setup_tissue();
void initialize_brain_delivery();
void advance_bbb_and_release(double time_min, double dt_min);
void begin_diffusion_mass_step();
void end_diffusion_mass_step();
void synchronize_endothelial_data();
void record_mass_balance(double time_min);
void write_mass_balance();

void brain_cell_phenotype(PhysiCell::Cell* cell, PhysiCell::Phenotype& phenotype, double dt);
std::vector<std::string> brain_delivery_coloring(PhysiCell::Cell* cell);

