#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"
#include "./custom_modules/brain_delivery.h"

using namespace BioFVM;
using namespace PhysiCell;

int main(int argc, char* argv[]) {
    try {
        const char* config = argc > 1 ? argv[1] : "./config/PhysiCell_settings.xml";
        if (!load_PhysiCell_config_file(config)) {
            return 2;
        }
        omp_set_num_threads(PhysiCell_settings.omp_num_threads);
        setup_microenvironment();
        create_cell_container_for_microenvironment(microenvironment, 30.0);
        create_cell_types();
        initialize_brain_delivery();
        setup_tissue();

        set_save_biofvm_mesh_as_matlab(true);
        set_save_biofvm_data_as_matlab(true);
        set_save_biofvm_cell_data(true);
        set_save_biofvm_cell_data_as_custom_matlab(true);

        char filename[1024];
        synchronize_endothelial_data();
        record_mass_balance(0.0);
        std::snprintf(filename, sizeof(filename), "%s/initial",
                      PhysiCell_settings.folder.c_str());
        save_PhysiCell_to_MultiCellDS_v2(filename, microenvironment, 0.0);

        BioFVM::RUNTIME_TIC();
        while (PhysiCell_globals.current_time <
               PhysiCell_settings.max_time + 0.1 * diffusion_dt) {
            if (PhysiCell_globals.current_time >
                PhysiCell_globals.next_full_save_time - 0.5 * diffusion_dt) {
                synchronize_endothelial_data();
                record_mass_balance(PhysiCell_globals.current_time);
                std::snprintf(filename, sizeof(filename), "%s/output%08u",
                              PhysiCell_settings.folder.c_str(),
                              PhysiCell_globals.full_output_index);
                save_PhysiCell_to_MultiCellDS_v2(
                    filename, microenvironment, PhysiCell_globals.current_time);
                ++PhysiCell_globals.full_output_index;
                PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
            }

            advance_bbb_and_release(PhysiCell_globals.current_time, diffusion_dt);
            begin_diffusion_mass_step();
            microenvironment.simulate_diffusion_decay(diffusion_dt);
            end_diffusion_mass_step();
            static_cast<Cell_Container*>(microenvironment.agent_container)
                ->update_all_cells(PhysiCell_globals.current_time);
            PhysiCell_globals.current_time += diffusion_dt;
        }

        synchronize_endothelial_data();
        record_mass_balance(PhysiCell_globals.current_time);
        std::snprintf(filename, sizeof(filename), "%s/final",
                      PhysiCell_settings.folder.c_str());
        save_PhysiCell_to_MultiCellDS_v2(
            filename, microenvironment, PhysiCell_globals.current_time);
        write_mass_balance();
        std::cout << "Completed brain-delivery simulation with " << all_cells->size()
                  << " cells.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "brain_delivery error: " << error.what() << '\n';
        return 2;
    }
}

