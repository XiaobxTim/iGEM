# PhysiCell–ParaView Brain Delivery Implementation Plan

> **For Codex:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a reproducible three-dose PhysiCell brain-delivery simulation and generate directly loadable ParaView 6.0.1 time-series data.

**Architecture:** Python exports Module 1–2 brain-blood boundary curves and orchestrates the runs. A Studio-compatible PhysiCell C++ project simulates BBB transfer, BioFVM AAV diffusion, fixed multicellular tissue, and per-cell Module 4–5 dynamics. A Python/VTK conversion stage creates validated VTI/VTP/VTM/PVD output.

**Tech Stack:** Python 3.12, NumPy, SciPy, PyYAML, pytest, PhysiCell 1.14.2, BioFVM, C++17/OpenMP, ParaView 6.0.1/pvpython.

---

### Task 1: PBPK boundary exporter

- Add failing tests for Module 1–2-only simulation, common normalization, CSV schema, and invalid input handling.
- Implement the exporter without modifying the existing full-model simulator.
- Verify the three dose curves and commit.

### Task 2: Deterministic tissue and parameter contract

- Add failing tests for the exact cell counts, cylindrical endothelial layout,
  non-overlapping parenchymal placement, rate conversion, and fixed seed.
- Implement geometry/config generation consumed by PhysiCell and ParaView.
- Verify reproducibility and commit.

### Task 3: Per-cell C++ Module 4–5

- Generate Python reference fixtures first.
- Implement the C++ state schema, fluxes, RK4 integration, positivity checks,
  and a standalone probe executable.
- Compare probe trajectories to Python with relative error below 1e-3 and
  commit.

### Task 4: PhysiCell/BioFVM project

- Add integration tests for configuration and expected raw output contract.
- Implement fixed agents, global BBB kinetics, perivascular source injection,
  cell uptake, intracellular callbacks, and mass accumulators.
- Build with pinned PhysiCell 1.14.2 and run short zero-dose and positive-dose
  smoke tests before committing.

### Task 5: ParaView conversion and metrics

- Add tests for metric definitions, VTK arrays, time-series manifests, and
  atomic output behavior.
- Implement compressed VTI/VTP/VTM/PVD writers using ParaView's pvpython and a
  three-dose comparison series.
- Round-trip all outputs through ParaView 6.0.1 and commit.

### Task 6: One-command pipeline and final verification

- Add the setup and orchestration commands, documentation, metadata, and
  overwrite protection.
- Run the full existing test suite, spatial tests, C++ build, three short
  scenarios, and ParaView round-trip checks.
- Run the 72-hour scenarios if performance permits; otherwise preserve the
  verified short run and report the exact external-runtime blocker.
