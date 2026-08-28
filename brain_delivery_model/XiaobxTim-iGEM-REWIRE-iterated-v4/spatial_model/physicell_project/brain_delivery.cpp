#include "brain_delivery.h"
#include "intracellular_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace BioFVM;
using namespace PhysiCell;

namespace {

struct BoundaryPoint {
    double time_min;
    double normalized;
};

struct BBBState {
    double surface = 0.0;
    double internalized = 0.0;
    double release_rate = 0.0;
};

struct BBBParameters {
    double blood_to_surface;
    double surface_to_endo;
    double surface_to_blood;
    double endo_to_isf;
    double endo_to_blood;
    double endo_loss;
};

struct MassRow {
    double time_min;
    double released;
    double extracellular;
    double cellular;
    double nuclear;
    double field_decay;
    double cell_loss;
    double nuclear_loss;
    double relative_error;
};

brain_delivery::Parameters intracellular_parameters{};
BBBParameters bbb_parameters{};
BBBState bbb_state{};
std::vector<BoundaryPoint> boundary_curve;
std::vector<int> release_voxels;
std::vector<MassRow> mass_rows;
int aav_index = -1;
double release_shell_volume = 0.0;
double source_scale = 1.0;
double cumulative_released = 0.0;
double cumulative_field_decay = 0.0;
double field_before_diffusion = 0.0;
double local_intracellular_dt = 1.0;

double numeric(const std::map<std::string, std::string>& values, const std::string& key) {
    const auto found = values.find(key);
    if (found == values.end()) {
        throw std::runtime_error("missing parameter: " + key);
    }
    const double value = std::stod(found->second);
    if (!std::isfinite(value) || value < 0.0) {
        throw std::runtime_error("invalid non-negative parameter: " + key);
    }
    return value;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) {
        result.push_back(value);
    }
    return result;
}

void load_boundary_curve(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open boundary CSV: " + path);
    }
    std::string line;
    if (!std::getline(input, line) ||
        line != "time_min,A_brain_blood_raw,A_brain_blood_normalized") {
        throw std::runtime_error("invalid boundary CSV header");
    }
    double previous = -1.0;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const auto values = split(line);
        if (values.size() != 3) {
            throw std::runtime_error("invalid boundary CSV row: " + line);
        }
        BoundaryPoint point{std::stod(values[0]), std::stod(values[2])};
        if (!std::isfinite(point.time_min) || !std::isfinite(point.normalized) ||
            point.time_min <= previous || point.normalized < 0.0) {
            throw std::runtime_error("boundary CSV must be finite, non-negative, and increasing");
        }
        boundary_curve.push_back(point);
        previous = point.time_min;
    }
    if (boundary_curve.size() < 2) {
        throw std::runtime_error("boundary CSV requires at least two rows");
    }
}

double boundary_value(double time_min) {
    if (time_min < boundary_curve.front().time_min || time_min > boundary_curve.back().time_min) {
        return 0.0;
    }
    const auto upper = std::lower_bound(
        boundary_curve.begin(), boundary_curve.end(), time_min,
        [](const BoundaryPoint& point, double time) { return point.time_min < time; });
    if (upper == boundary_curve.begin()) {
        return upper->normalized;
    }
    if (upper == boundary_curve.end()) {
        return boundary_curve.back().normalized;
    }
    const auto lower = upper - 1;
    const double fraction = (time_min - lower->time_min) / (upper->time_min - lower->time_min);
    return lower->normalized + fraction * (upper->normalized - lower->normalized);
}

std::array<double, 2> bbb_rhs(const std::array<double, 2>& state, double blood) {
    return {
        bbb_parameters.blood_to_surface * blood -
            (bbb_parameters.surface_to_endo + bbb_parameters.surface_to_blood) * state[0],
        bbb_parameters.surface_to_endo * state[0] -
            (bbb_parameters.endo_to_isf + bbb_parameters.endo_to_blood +
             bbb_parameters.endo_loss) * state[1],
    };
}

void update_bbb(double time_min, double dt_min) {
    const std::array<double, 2> state{bbb_state.surface, bbb_state.internalized};
    const auto k1 = bbb_rhs(state, boundary_value(time_min));
    const std::array<double, 2> state2{state[0] + 0.5 * dt_min * k1[0],
                                       state[1] + 0.5 * dt_min * k1[1]};
    const auto k2 = bbb_rhs(state2, boundary_value(time_min + 0.5 * dt_min));
    const std::array<double, 2> state3{state[0] + 0.5 * dt_min * k2[0],
                                       state[1] + 0.5 * dt_min * k2[1]};
    const auto k3 = bbb_rhs(state3, boundary_value(time_min + 0.5 * dt_min));
    const std::array<double, 2> state4{state[0] + dt_min * k3[0],
                                       state[1] + dt_min * k3[1]};
    const auto k4 = bbb_rhs(state4, boundary_value(time_min + dt_min));
    bbb_state.surface = std::max(
        state[0] + dt_min * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0, 0.0);
    bbb_state.internalized = std::max(
        state[1] + dt_min * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0, 0.0);
    bbb_state.release_rate = bbb_parameters.endo_to_isf * bbb_state.internalized;
}

double field_mass() {
    double total = 0.0;
    const int voxel_count = static_cast<int>(microenvironment.number_of_voxels());
    for (int index = 0; index < voxel_count; ++index) {
        total += microenvironment.density_vector(index)[aav_index] *
                 microenvironment.mesh.voxels[index].volume;
    }
    return total;
}

double custom_sum(const std::string& name) {
    double total = 0.0;
    for (Cell* cell : *all_cells) {
        if (cell->type == 2 || cell->type == 3) {
            total += cell->custom_data[name];
        }
    }
    return total;
}

void initialize_cell_state(Cell* cell, double apoe_scale, double distance) {
    const auto state = brain_delivery::initial_state(intracellular_parameters, apoe_scale);
    for (std::size_t index = 0; index < state.size(); ++index) {
        cell->custom_data[brain_delivery::state_names()[index]] = state[index];
    }
    cell->custom_data["distance_to_vessel_um"] = distance;
    cell->custom_data["editing_fraction"] = 0.0;
    cell->custom_data["off_target_burden"] = 0.0;
    cell->custom_data["vector_cell_loss_cumulative"] = 0.0;
    cell->custom_data["vector_nuclear_loss_cumulative"] = 0.0;
}

}  // namespace

void setup_microenvironment() {
    initialize_microenvironment();
}

void create_cell_types() {
    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);
    cell_defaults.functions.volume_update_function = nullptr;
    cell_defaults.functions.update_velocity = nullptr;
    cell_defaults.functions.update_migration_bias = nullptr;
    cell_defaults.functions.update_phenotype = nullptr;
    cell_defaults.functions.custom_cell_rule = nullptr;
    cell_defaults.functions.contact_function = nullptr;
    initialize_cell_definitions_from_pugixml();
    build_cell_definitions_maps();
    setup_signal_behavior_dictionaries();
    for (const std::string name : {"default", "endothelial", "neuron", "astrocyte"}) {
        Cell_Definition* definition = find_cell_definition(name);
        definition->functions.volume_update_function = nullptr;
        definition->functions.update_velocity = nullptr;
        definition->functions.custom_cell_rule = nullptr;
        definition->functions.contact_function = nullptr;
        definition->functions.update_phenotype =
            (name == "neuron" || name == "astrocyte") ? brain_cell_phenotype : nullptr;
    }
}

void setup_tissue() {
    const std::string path = parameters.strings("cells_csv");
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open cell CSV: " + path);
    }
    std::string line;
    std::getline(input, line);
    int count = 0;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        const auto values = split(line);
        if (values.size() != 7) {
            throw std::runtime_error("invalid cell CSV row: " + line);
        }
        const std::string type = values[1];
        Cell_Definition* definition = find_cell_definition(type);
        if (definition == nullptr) {
            throw std::runtime_error("unknown cell type in CSV: " + type);
        }
        Cell* cell = create_cell(*definition);
        const double radius = std::stod(values[5]);
        cell->set_total_volume(4.0 * M_PI * radius * radius * radius / 3.0);
        cell->assign_position({std::stod(values[2]), std::stod(values[3]), std::stod(values[4])});
        cell->is_movable = false;
        const double distance = std::stod(values[6]);
        if (type == "neuron") {
            initialize_cell_state(cell, 0.2, distance);
        } else if (type == "astrocyte") {
            initialize_cell_state(cell, 1.0, distance);
        } else {
            cell->custom_data["distance_to_vessel_um"] = distance;
        }
        ++count;
    }
    if (count != 1480 && parameters.ints("random_seed") == 42) {
        std::cout << "Loaded " << count << " cells (non-default geometry count)." << std::endl;
    }
}

void initialize_brain_delivery() {
    aav_index = microenvironment.find_density_index("extracellular_AAV");
    if (aav_index < 0) {
        throw std::runtime_error("extracellular_AAV substrate is missing");
    }
    intracellular_parameters = brain_delivery::load_parameters(parameters.strings("parameter_file"));
    const auto flat = brain_delivery::read_parameter_file(parameters.strings("parameter_file"));
    bbb_parameters = {
        numeric(flat, "bbb.k_brainblood_to_EC_per_min"),
        numeric(flat, "bbb.k_EC_to_endo_per_min"),
        numeric(flat, "bbb.k_EC_to_brainblood_per_min"),
        numeric(flat, "bbb.k_endo_to_ISF_per_min"),
        numeric(flat, "bbb.k_endo_to_brainblood_per_min"),
        numeric(flat, "bbb.k_endo_loss_per_min"),
    };
    source_scale = parameters.doubles("source_scale");
    local_intracellular_dt = parameters.doubles("intracellular_dt_min");
    load_boundary_curve(parameters.strings("boundary_csv"));

    const double inner = parameters.doubles("vessel_radius_um") +
                         parameters.doubles("endothelial_radius_um");
    const double outer = inner + parameters.doubles("perivascular_shell_thickness_um");
    const int voxel_count = static_cast<int>(microenvironment.number_of_voxels());
    for (int index = 0; index < voxel_count; ++index) {
        const auto& center = microenvironment.mesh.voxels[index].center;
        const double radial = std::hypot(center[1], center[2]);
        if (radial >= inner && radial < outer) {
            release_voxels.push_back(index);
            release_shell_volume += microenvironment.mesh.voxels[index].volume;
        }
    }
    if (release_voxels.empty() || release_shell_volume <= 0.0) {
        throw std::runtime_error("perivascular release shell contains no voxels");
    }
}

void advance_bbb_and_release(double time_min, double dt_min) {
    update_bbb(time_min, dt_min);
    const double amount = bbb_state.release_rate * source_scale * dt_min;
    const double concentration_increment = amount / release_shell_volume;
    for (int voxel : release_voxels) {
        microenvironment.density_vector(voxel)[aav_index] += concentration_increment;
    }
    cumulative_released += amount;
}

void begin_diffusion_mass_step() {
    field_before_diffusion = field_mass();
}

void end_diffusion_mass_step() {
    cumulative_field_decay += std::max(field_before_diffusion - field_mass(), 0.0);
}

void brain_cell_phenotype(Cell* cell, Phenotype& phenotype, double dt) {
    brain_delivery::State state{};
    state[brain_delivery::A_CELL] = phenotype.molecular.internalized_total_substrates[aav_index];
    for (std::size_t index = 1; index < state.size(); ++index) {
        state[index] = cell->custom_data[brain_delivery::state_names()[index]];
    }
    const double apoe_scale = cell->type == 3 ? 1.0 : 0.2;
    const int substeps = std::max(
        1, static_cast<int>(std::ceil(dt / local_intracellular_dt)));
    const double sub_dt = dt / substeps;
    for (int step = 0; step < substeps; ++step) {
        const auto previous = state;
        state = brain_delivery::rk4_step(
            state, intracellular_parameters, 0.0, apoe_scale, sub_dt);
        cell->custom_data["vector_cell_loss_cumulative"] +=
            0.5 * intracellular_parameters.k_cell_loss *
            (previous[brain_delivery::A_CELL] + state[brain_delivery::A_CELL]) * sub_dt;
        cell->custom_data["vector_nuclear_loss_cumulative"] +=
            0.5 * intracellular_parameters.k_deg_v *
            (previous[brain_delivery::A_NUCLEUS] + state[brain_delivery::A_NUCLEUS]) * sub_dt;
    }
    phenotype.molecular.internalized_total_substrates[aav_index] = state[brain_delivery::A_CELL];
    for (std::size_t index = 0; index < state.size(); ++index) {
        cell->custom_data[brain_delivery::state_names()[index]] = state[index];
    }
    const double total_apoe = state[brain_delivery::S_APOE4] +
        state[brain_delivery::S_APOE3_LIKE] + state[brain_delivery::S_APOE2_LIKE] +
        state[brain_delivery::S_APOE158_ONLY] + state[brain_delivery::C_APOE112] +
        state[brain_delivery::C_APOE158] + state[brain_delivery::C_APOE158_AFTER112] +
        state[brain_delivery::C_APOE112_AFTER158] + 1e-12;
    const double therapeutic =
        state[brain_delivery::S_APOE3_LIKE] + state[brain_delivery::S_APOE2_LIKE];
    const double off_target = state[brain_delivery::B_LOCAL_BYSTANDER] +
        state[brain_delivery::E_PUF_OFF] + state[brain_delivery::E_DEAMINASE_BG];
    cell->custom_data["editing_fraction"] = therapeutic / total_apoe;
    cell->custom_data["off_target_burden"] = off_target / total_apoe;
    cell->custom_data["uptake_rate"] = phenotype.secretion.uptake_rates[aav_index];
}

void synchronize_endothelial_data() {
    for (Cell* cell : *all_cells) {
        if (cell->type == 1) {
            cell->custom_data["endothelial_surface_AAV"] = bbb_state.surface;
            cell->custom_data["endothelial_internalized_AAV"] = bbb_state.internalized;
            cell->custom_data["BBB_release_rate"] = bbb_state.release_rate;
        }
    }
}

void record_mass_balance(double time_min) {
    double cellular = 0.0;
    double nuclear = 0.0;
    for (Cell* cell : *all_cells) {
        if (cell->type == 2 || cell->type == 3) {
            cellular += cell->phenotype.molecular.internalized_total_substrates[aav_index];
            nuclear += cell->custom_data["A_nucleus"];
        }
    }
    const double cell_loss = custom_sum("vector_cell_loss_cumulative");
    const double nuclear_loss = custom_sum("vector_nuclear_loss_cumulative");
    const double extracellular = field_mass();
    const double accounted = extracellular + cellular + nuclear + cumulative_field_decay +
                             cell_loss + nuclear_loss;
    const double relative_error = cumulative_released > 1e-15
        ? std::abs(cumulative_released - accounted) / cumulative_released
        : std::abs(accounted);
    mass_rows.push_back({time_min, cumulative_released, extracellular, cellular, nuclear,
                         cumulative_field_decay, cell_loss, nuclear_loss, relative_error});
}

void write_mass_balance() {
    const std::string path = PhysiCell_settings.folder + "/mass_balance.csv";
    std::ofstream output(path);
    output << "time_min,released,extracellular,cellular,nuclear,field_decay,cell_loss,nuclear_loss,relative_error\n";
    output << std::setprecision(17);
    for (const auto& row : mass_rows) {
        output << row.time_min << ',' << row.released << ',' << row.extracellular << ','
               << row.cellular << ',' << row.nuclear << ',' << row.field_decay << ','
               << row.cell_loss << ',' << row.nuclear_loss << ',' << row.relative_error << '\n';
    }
}

std::vector<std::string> brain_delivery_coloring(Cell* cell) {
    if (cell->type == 1) return {"cyan", "cyan", "cyan", "black"};
    if (cell->type == 2) return {"royalblue", "royalblue", "royalblue", "black"};
    if (cell->type == 3) return {"orange", "orange", "orange", "black"};
    return {"grey", "grey", "grey", "black"};
}
