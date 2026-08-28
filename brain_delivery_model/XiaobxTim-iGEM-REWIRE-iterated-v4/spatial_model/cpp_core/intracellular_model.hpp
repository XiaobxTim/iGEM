#pragma once

#include <array>
#include <map>
#include <string>
#include <vector>

namespace brain_delivery {

enum StateIndex {
    A_CELL,
    A_NUCLEUS,
    EDITOR_MRNA,
    EDITOR_PROTEIN,
    S_APOE4,
    S_APOE3_LIKE,
    S_APOE2_LIKE,
    S_APOE158_ONLY,
    C_APOE112,
    C_APOE158,
    C_APOE158_AFTER112,
    C_APOE112_AFTER158,
    B_LOCAL_BYSTANDER,
    S_PUF_OFF,
    C_PUF_OFF,
    E_PUF_OFF,
    S_DEAMINASE_BG,
    E_DEAMINASE_BG,
    STATE_COUNT
};

using State = std::array<double, STATE_COUNT>;

struct Parameters {
    double k_cell_to_nuc;
    double k_cell_loss;
    double k_deg_v;
    double k_tx;
    double k_deg_m;
    double k_tl;
    double k_deg_p;
    double s_apoe4_init;
    double s_puf_off_init;
    double s_deaminase_bg_init;
    double k_on_112;
    double k_off_112;
    double k_cat_112;
    double k_on_158;
    double k_off_158;
    double k_cat_158;
    double k_prod_apoe;
    double k_deg_apoe;
    double local_bystander_per_112;
    double local_bystander_per_158;
    double k_on_puf_off;
    double k_off_puf_off;
    double k_cat_puf_off;
    double k_deaminase_bg;
    double k_prod_puf_off;
    double k_deg_puf_off;
    double k_prod_deaminase_bg;
    double k_deg_deaminase_bg;
    double cat_scale;
    double uc_context_scale;
    double background_scale;
};

std::map<std::string, std::string> read_parameter_file(const std::string& path);
Parameters load_parameters(const std::string& path);
const std::array<const char*, STATE_COUNT>& state_names();
State initial_state(const Parameters& parameters, double apoe_scale);
State rhs(const State& state, const Parameters& parameters, double uptake, double apoe_scale);
State rk4_step(const State& state, const Parameters& parameters, double uptake,
               double apoe_scale, double dt_min);

}  // namespace brain_delivery

