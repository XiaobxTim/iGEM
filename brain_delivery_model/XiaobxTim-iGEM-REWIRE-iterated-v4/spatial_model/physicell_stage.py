from __future__ import annotations

import shutil
from pathlib import Path


def _render_makefile(template: str) -> str:
    object_marker = "BioFVM_OBJECTS :="
    marker_index = template.index(object_marker)
    prefix = """VERSION := $(shell grep . VERSION.txt | cut -f1 -d:)
PROGRAM_NAME := brain_delivery

CC := clang++
ifdef PHYSICELL_CPP
    CC := $(PHYSICELL_CPP)
endif

CFLAGS := -O3 -fomit-frame-pointer -m64 -std=c++17 -Wall -Wextra
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CFLAGS += -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
    LDLIBS := -L/opt/homebrew/opt/libomp/lib -lomp
else
    CFLAGS += -fopenmp
    LDLIBS := -fopenmp
endif

CFLAGS_LINK := $(CFLAGS)
COMPILE_COMMAND := $(CC) $(CFLAGS) $(EXTRA_FLAGS)
LINK_COMMAND := $(CC) $(CFLAGS_LINK) $(EXTRA_FLAGS)

"""
    rendered = prefix + template[marker_index:]
    rendered = rendered.replace(
        "PhysiCell_custom_module_OBJECTS := custom.o",
        "PhysiCell_custom_module_OBJECTS := brain_delivery.o intracellular_model.o",
    )
    rendered = rendered.replace(
        "$(COMPILE_COMMAND) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp",
        "$(COMPILE_COMMAND) -o $(PROGRAM_NAME) $(ALL_OBJECTS) main.cpp $(LDLIBS)",
    )
    rendered = rendered.replace(
        "custom.o: ./custom_modules/custom.cpp \n\t$(COMPILE_COMMAND) -c ./custom_modules/custom.cpp",
        "brain_delivery.o: ./custom_modules/brain_delivery.cpp ./custom_modules/brain_delivery.h\n"
        "\t$(COMPILE_COMMAND) -c ./custom_modules/brain_delivery.cpp\n\n"
        "intracellular_model.o: ./custom_modules/intracellular_model.cpp ./custom_modules/intracellular_model.hpp\n"
        "\t$(COMPILE_COMMAND) -c ./custom_modules/intracellular_model.cpp",
    )
    return rendered


def stage_physicell_project(physicell_root: str | Path, spatial_model_root: str | Path) -> Path:
    """Install the tracked project overlay into a staged PhysiCell 1.14.2 tree."""

    root = Path(physicell_root)
    spatial = Path(spatial_model_root)
    version = root / "VERSION.txt"
    if not version.exists() or version.read_text(encoding="utf-8").strip() != "1.14.2":
        raise ValueError(f"expected a PhysiCell 1.14.2 tree at {root}")
    sample_makefile = root / "sample_projects" / "template" / "Makefile"
    if not sample_makefile.exists():
        raise FileNotFoundError(sample_makefile)

    project = spatial / "physicell_project"
    custom_destination = root / "custom_modules"
    custom_destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(project / "main.cpp", root / "main.cpp")
    shutil.copy2(project / "brain_delivery.cpp", custom_destination / "brain_delivery.cpp")
    shutil.copy2(project / "brain_delivery.h", custom_destination / "brain_delivery.h")
    shutil.copy2(spatial / "cpp_core" / "intracellular_model.cpp", custom_destination / "intracellular_model.cpp")
    shutil.copy2(spatial / "cpp_core" / "intracellular_model.hpp", custom_destination / "intracellular_model.hpp")
    (root / "Makefile").write_text(
        _render_makefile(sample_makefile.read_text(encoding="utf-8")),
        encoding="utf-8",
    )
    return root

