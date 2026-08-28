CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic

.PHONY: all clean

all: ode_probe

ode_probe: ode_probe.cpp intracellular_model.cpp intracellular_model.hpp
	$(CXX) $(CXXFLAGS) ode_probe.cpp intracellular_model.cpp -o $@

clean:
	rm -f ode_probe
