#include "intracellular_model.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::string argument(int argc, char** argv, const std::string& name,
                     const std::string& fallback = "") {
    for (int index = 1; index + 1 < argc; ++index) {
        if (argv[index] == name) {
            return argv[index + 1];
        }
    }
    return fallback;
}

double numeric_argument(int argc, char** argv, const std::string& name,
                        double fallback) {
    const auto value = argument(argc, argv, name);
    return value.empty() ? fallback : std::stod(value);
}

void print_row(double time, const brain_delivery::State& state) {
    std::cout << std::setprecision(17) << time;
    for (double value : state) {
        std::cout << ',' << value;
    }
    std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const auto parameter_path = argument(argc, argv, "--parameters");
        if (parameter_path.empty()) {
            throw std::runtime_error("--parameters is required");
        }
        const double minutes = numeric_argument(argc, argv, "--minutes", 120.0);
        const double dt = numeric_argument(argc, argv, "--dt", 1.0);
        const double uptake = numeric_argument(argc, argv, "--uptake", 0.0005);
        const double apoe_scale = numeric_argument(argc, argv, "--apoe-scale", 1.0);
        if (!std::isfinite(uptake) || uptake < 0.0) {
            throw std::runtime_error("uptake must be non-negative");
        }
        if (!std::isfinite(minutes) || minutes <= 0.0 || !std::isfinite(dt) ||
            dt <= 0.0 || std::abs(std::round(minutes / dt) - minutes / dt) > 1e-9) {
            throw std::runtime_error("minutes must be a positive multiple of dt");
        }

        const auto parameters = brain_delivery::load_parameters(parameter_path);
        auto state = brain_delivery::initial_state(parameters, apoe_scale);
        std::cout << "time_min";
        for (const auto* name : brain_delivery::state_names()) {
            std::cout << ',' << name;
        }
        std::cout << '\n';
        print_row(0.0, state);
        const int steps = static_cast<int>(std::llround(minutes / dt));
        for (int step = 1; step <= steps; ++step) {
            state = brain_delivery::rk4_step(state, parameters, uptake, apoe_scale, dt);
            print_row(step * dt, state);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}

