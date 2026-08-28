# PhysiCell–ParaView Brain Delivery Design

## Goal

Extend the existing well-mixed AAV delivery model with a three-dimensional,
agent-based brain microvascular tissue model for iGEM demonstration. The model
must preserve the existing Module 1–5 mechanisms while adding BBB release,
interstitial diffusion, cell-specific uptake, and spatially heterogeneous APOE
editing.

## Architecture

The existing Python Module 1–2 equations generate a prescribed brain-blood AAV
time course for doses 0.3, 1.0, and 3.0. PhysiCell 1.14.2 reads those curves,
integrates the existing BBB kinetics, releases normalized AAV into a BioFVM
field around a straight capillary, and integrates full Module 4–5 kinetics in
each neuron and astrocyte. The coupling is one-way: the spatial simulation does
not feed back into the PBPK compartments.

The tissue domain is 400 x 300 x 300 micrometres with a 10 micrometre mesh. A
12 micrometre radius vessel runs along the x-axis and is represented by 280
fixed endothelial agents. The parenchyma contains 800 neurons and 400
astrocytes. All three dose scenarios use identical cell locations generated
with seed 42.

The sole diffusing substrate is normalized extracellular AAV. It diffuses at
10 square micrometres per minute, decays with a 24-hour half-life, is released
in a 10 micrometre perivascular shell, and is taken up by parenchymal agents.
Cells are non-motile, non-dividing, and non-dying in this first version.

## Intracellular model

Each neuron and astrocyte tracks cellular and nuclear AAV, editor mRNA and
protein, APOE4/APOE3-like/APOE2-like/mixed-edit substrates and products, four
enzyme-substrate complexes, and all three existing off-target classes. Existing
hour-based parameters are converted to minute-based rates. The C++ model uses
RK4 with one-minute substeps and must match the Python reference trajectory to
relative error below 1e-3.

Astrocytes use 1.4 times the neuronal AAV uptake and the existing APOE substrate
and production scale. Neurons use 0.2 times the astrocyte APOE scale. These are
configurable demonstration priors, not calibrated in-vivo estimates.

## Data products

Each dose produces raw MultiCellDS snapshots, metrics, mass balance, metadata,
and a ParaView time series. Structured fields are written as VTI, cell points
and attributes as VTP, each frame as VTM, and the time series as PVD. A combined
three-dose PVD places the scenarios side-by-side. ParaView 6.0.1's bundled
`pvpython` writes and validates the compressed VTK XML files.

## Validation and interpretation

Zero-dose and zero-BBB controls must stay at zero. Released vector mass must be
balanced among extracellular, cellular, nuclear, degraded, and lost pools to
within one percent. Tissue AAV exposure and average editing must increase from
low to medium to high dose. Near-vessel cells must exceed far-vessel cells on
average, and astrocytes must exceed neurons under the configured priors.

Geometry and time have physical units; AAV abundance remains normalized. The
model is a reaction-diffusion and multicellular demonstration, not a fluid-flow
or calibrated clinical prediction model.
