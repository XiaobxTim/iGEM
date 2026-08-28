#include "intracellular_model.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace brain_delivery {
namespace {

double number(const std::map<std::string, std::string>& values, const std::string& key) {
    auto found = values.find(key);
    if (found == values.end()) {
        throw std::runtime_error("missing parameter: " + key);
    }
    std::size_t parsed = 0;
    const double value = std::stod(found->second, &parsed);
    if (parsed != found->second.size() || !std::isfinite(value)) {
        throw std::runtime_error("invalid numeric parameter: " + key);
    }
    return value;
}

std::string text_value(const std::map<std::string, std::string>& values,
                       const std::string& key) {
    auto found = values.find(key);
    if (found == values.end()) {
        throw std::runtime_error("missing parameter: " + key);
    }
    return found->second;
}

State add_scaled(const State& left, const State& right, double scale) {
    State result{};
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] = left[index] + scale * right[index];
    }
    return result;
}

}  // namespace

std::map<std::string, std::string> read_parameter_file(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open parameter file: " + path);
    }
    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') {
            continue;
        }
        const auto separator = line.find('=');
        if (separator == std::string::npos || separator == 0) {
            throw std::runtime_error("invalid parameter line: " + line);
        }
        values[line.substr(0, separator)] = line.substr(separator + 1);
    }
    return values;
}

Parameters load_parameters(const std::string& path) {
    const auto values = read_parameter_file(path);
    Parameters p{};
    p.k_cell_to_nuc = number(values, "intracellular.k_cell_to_nuc_per_min");
    p.k_cell_loss = number(values, "intracellular.k_cell_loss_per_min");
    p.k_deg_v = number(values, "intracellular.k_deg_v_per_min");
    p.k_tx = number(values, "intracellular.k_tx_per_min");
    p.k_deg_m = number(values, "intracellular.k_deg_m_per_min");
    p.k_tl = number(values, "intracellular.k_tl_per_min");
    p.k_deg_p = number(values, "intracellular.k_deg_p_per_min");
    p.s_apoe4_init = number(values, "editing.S_APOE4_init");
    p.s_puf_off_init = number(values, "editing.S_puf_off_init");
    p.s_deaminase_bg_init = number(values, "editing.S_deaminase_bg_init");
    p.k_on_112 = number(values, "editing.k_on_112_per_min");
    p.k_off_112 = number(values, "editing.k_off_112_per_min");
    p.k_cat_112 = number(values, "editing.k_cat_112_per_min");
    p.k_on_158 = number(values, "editing.k_on_158_per_min");
    p.k_off_158 = number(values, "editing.k_off_158_per_min");
    p.k_cat_158 = number(values, "editing.k_cat_158_per_min");
    p.k_prod_apoe = number(values, "editing.k_prod_apoe_per_min");
    p.k_deg_apoe = number(values, "editing.k_deg_apoe_per_min");
    p.local_bystander_per_112 = number(values, "editing.local_bystander_per_112");
    p.local_bystander_per_158 = number(values, "editing.local_bystander_per_158");
    p.k_on_puf_off = number(values, "editing.k_on_puf_off_per_min");
    p.k_off_puf_off = number(values, "editing.k_off_puf_off_per_min");
    p.k_cat_puf_off = number(values, "editing.k_cat_puf_off_per_min");
    p.k_deaminase_bg = number(values, "editing.k_deaminase_bg_per_min");
    p.k_prod_puf_off = number(values, "editing.k_prod_puf_off_per_min");
    p.k_deg_puf_off = number(values, "editing.k_deg_puf_off_per_min");
    p.k_prod_deaminase_bg = number(values, "editing.k_prod_deaminase_bg_per_min");
    p.k_deg_deaminase_bg = number(values, "editing.k_deg_deaminase_bg_per_min");

    const auto editor = text_value(values, "editing.editor_type");
    if (editor == "A3A" || editor == "a3a" || editor == "APOBEC3A") {
        p.cat_scale = 1.0;
        p.uc_context_scale = number(values, "editing.uc_context_scale");
        p.background_scale = number(values, "editing.a3a_background_scale");
    } else if (editor == "APOBEC1" || editor == "apobec1") {
        p.cat_scale = number(values, "editing.apobec1_cat_scale");
        p.uc_context_scale = 0.9;
        p.background_scale = number(values, "editing.apobec1_background_scale");
    } else if (editor == "ProAPOBEC" || editor == "proapobec") {
        p.cat_scale = number(values, "editing.proapobec_cat_scale");
        p.uc_context_scale = 1.1;
        p.background_scale = number(values, "editing.proapobec_background_scale");
    } else {
        throw std::runtime_error("unsupported editor_type: " + editor);
    }
    return p;
}

const std::array<const char*, STATE_COUNT>& state_names() {
    static const std::array<const char*, STATE_COUNT> names = {
        "A_cell", "A_nucleus", "editor_mRNA", "editor_protein", "S_APOE4",
        "S_APOE3_like", "S_APOE2_like", "S_APOE158_only", "C_APOE112",
        "C_APOE158", "C_APOE158_after112", "C_APOE112_after158",
        "B_local_bystander", "S_puf_off", "C_puf_off", "E_puf_off",
        "S_deaminase_bg", "E_deaminase_bg"};
    return names;
}

State initial_state(const Parameters& p, double apoe_scale) {
    if (!std::isfinite(apoe_scale) || apoe_scale < 0.0) {
        throw std::runtime_error("apoe_scale must be non-negative");
    }
    State state{};
    state[S_APOE4] = p.s_apoe4_init * apoe_scale;
    state[S_PUF_OFF] = p.s_puf_off_init;
    state[S_DEAMINASE_BG] = p.s_deaminase_bg_init;
    return state;
}

State rhs(const State& s, const Parameters& p, double uptake, double apoe_scale) {
    if (!std::isfinite(uptake) || uptake < 0.0) {
        throw std::runtime_error("uptake must be non-negative");
    }
    State d{};
    d[A_CELL] = uptake - (p.k_cell_to_nuc + p.k_cell_loss) * s[A_CELL];
    d[A_NUCLEUS] = p.k_cell_to_nuc * s[A_CELL] - p.k_deg_v * s[A_NUCLEUS];
    d[EDITOR_MRNA] = p.k_tx * s[A_NUCLEUS] - p.k_deg_m * s[EDITOR_MRNA];
    d[EDITOR_PROTEIN] = p.k_tl * s[EDITOR_MRNA] - p.k_deg_p * s[EDITOR_PROTEIN];

    const double protein = std::max(s[EDITOR_PROTEIN], 0.0);
    const double cat112 = p.k_cat_112 * p.cat_scale * p.uc_context_scale;
    const double cat158 = p.k_cat_158 * p.cat_scale;
    const double cat_off = p.k_cat_puf_off * p.cat_scale;
    const double k_bg = p.k_deaminase_bg * p.background_scale;

    const double bind112 = p.k_on_112 * protein * std::max(s[S_APOE4], 0.0);
    const double unbind112 = p.k_off_112 * std::max(s[C_APOE112], 0.0);
    const double edit112 = cat112 * std::max(s[C_APOE112], 0.0);
    const double bind158 = p.k_on_158 * protein * std::max(s[S_APOE4], 0.0);
    const double unbind158 = p.k_off_158 * std::max(s[C_APOE158], 0.0);
    const double edit158 = cat158 * std::max(s[C_APOE158], 0.0);
    const double bind158_after112 =
        p.k_on_158 * protein * std::max(s[S_APOE3_LIKE], 0.0);
    const double unbind158_after112 =
        p.k_off_158 * std::max(s[C_APOE158_AFTER112], 0.0);
    const double edit158_after112 =
        cat158 * std::max(s[C_APOE158_AFTER112], 0.0);
    const double bind112_after158 =
        p.k_on_112 * protein * std::max(s[S_APOE158_ONLY], 0.0);
    const double unbind112_after158 =
        p.k_off_112 * std::max(s[C_APOE112_AFTER158], 0.0);
    const double edit112_after158 =
        cat112 * std::max(s[C_APOE112_AFTER158], 0.0);
    const double bind_puf_off =
        p.k_on_puf_off * protein * std::max(s[S_PUF_OFF], 0.0);
    const double unbind_puf_off = p.k_off_puf_off * std::max(s[C_PUF_OFF], 0.0);
    const double edit_puf_off = cat_off * std::max(s[C_PUF_OFF], 0.0);
    const double deaminase_bg_edit = k_bg * protein * std::max(s[S_DEAMINASE_BG], 0.0);
    const double local_bystander =
        p.local_bystander_per_112 * (edit112 + edit112_after158) +
        p.local_bystander_per_158 * (edit158 + edit158_after112);

    const double binding_loss = bind112 + bind158 + bind158_after112 +
                                bind112_after158 + bind_puf_off;
    const double editor_release = unbind112 + unbind158 + unbind158_after112 +
                                  unbind112_after158 + unbind_puf_off + edit112 +
                                  edit158 + edit158_after112 + edit112_after158 +
                                  edit_puf_off;
    d[EDITOR_PROTEIN] += -binding_loss + editor_release;

    d[S_APOE4] = p.k_prod_apoe * apoe_scale - p.k_deg_apoe * s[S_APOE4] -
                 bind112 + unbind112 - bind158 + unbind158;
    d[C_APOE112] = bind112 - unbind112 - edit112;
    d[C_APOE158] = bind158 - unbind158 - edit158;
    d[S_APOE3_LIKE] = edit112 - p.k_deg_apoe * s[S_APOE3_LIKE] -
                      bind158_after112 + unbind158_after112;
    d[C_APOE158_AFTER112] = bind158_after112 - unbind158_after112 - edit158_after112;
    d[S_APOE158_ONLY] = edit158 - p.k_deg_apoe * s[S_APOE158_ONLY] -
                        bind112_after158 + unbind112_after158;
    d[C_APOE112_AFTER158] = bind112_after158 - unbind112_after158 - edit112_after158;
    d[S_APOE2_LIKE] = edit158_after112 + edit112_after158 -
                      p.k_deg_apoe * s[S_APOE2_LIKE];
    d[B_LOCAL_BYSTANDER] = local_bystander;
    d[S_PUF_OFF] = p.k_prod_puf_off - p.k_deg_puf_off * s[S_PUF_OFF] -
                   bind_puf_off + unbind_puf_off;
    d[C_PUF_OFF] = bind_puf_off - unbind_puf_off - edit_puf_off;
    d[E_PUF_OFF] = edit_puf_off;
    d[S_DEAMINASE_BG] = p.k_prod_deaminase_bg -
                        p.k_deg_deaminase_bg * s[S_DEAMINASE_BG] -
                        deaminase_bg_edit;
    d[E_DEAMINASE_BG] = deaminase_bg_edit;
    return d;
}

State rk4_step(const State& state, const Parameters& parameters, double uptake,
               double apoe_scale, double dt_min) {
    if (!std::isfinite(dt_min) || dt_min <= 0.0) {
        throw std::runtime_error("dt must be positive");
    }
    const State k1 = rhs(state, parameters, uptake, apoe_scale);
    const State k2 = rhs(add_scaled(state, k1, 0.5 * dt_min), parameters, uptake,
                         apoe_scale);
    const State k3 = rhs(add_scaled(state, k2, 0.5 * dt_min), parameters, uptake,
                         apoe_scale);
    const State k4 = rhs(add_scaled(state, k3, dt_min), parameters, uptake, apoe_scale);
    State updated{};
    for (std::size_t index = 0; index < updated.size(); ++index) {
        updated[index] = state[index] + dt_min *
            (k1[index] + 2.0 * k2[index] + 2.0 * k3[index] + k4[index]) / 6.0;
        if (!std::isfinite(updated[index]) || updated[index] < -1e-9) {
            throw std::runtime_error(std::string("invalid intracellular state: ") +
                                     state_names()[index]);
        }
        updated[index] = std::max(updated[index], 0.0);
    }
    return updated;
}

}  // namespace brain_delivery

