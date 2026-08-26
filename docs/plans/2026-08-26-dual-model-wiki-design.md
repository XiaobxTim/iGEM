# Dual-Model iGEM Wiki Design

## Purpose

Create an English, model-focused iGEM Wiki that presents the Brain Delivery
Digital Twin v4 and PUF-OffTarget Atlas as two complementary but independent
models. The Wiki is static; each model has a separately deployable interactive
application.

## Information architecture

The Wiki contains Home, Model, Brain Delivery, Off-Target Atlas, Engineering,
Software, and Resources pages. Its core visual is a dual-track systems map:
AAV delivery runs from administration to brain expression while PUF design runs
from recognition sequence to transcriptome-wide candidate prioritization. The
tracks meet at safer RNA-editing design.

## Visual direction

Use a scientific field-atlas aesthetic: navy structure, ivory graph-paper
surfaces, editorial typography, cobalt/amber for Brain Delivery, teal/coral for
PUF Atlas, and the established A/C/G/U palette. Fonts and assets are self-hosted.
The experience is responsive, keyboard accessible, and reduced-motion aware.

## Application architecture

The Wiki uses the iGEM React/Vite pattern and configurable team/App URLs. PUF
Atlas remains a FastAPI application. Brain Delivery gains a lightweight
FastAPI/Jinja/Plotly interface over the existing v4 simulator. The applications
exchange data through a downloadable/uploadable candidate-panel CSV rather than
server-to-server calls.

## Scientific guardrails

Content is derived from checked-in v4 Wiki/DBTL/parameter documents and PUF
documentation. Predictions are consistently labeled simulated or prioritized,
not experimentally proven or clinically calibrated. Missing PUF expression or
accessibility evidence maps conservatively to 1.0 in the Brain panel and is
recorded in metadata.

## Acceptance

Seven English Wiki pages must build without external runtime assets. Both Apps
must run locally. A synthetic PUF scan must export a candidate panel that the
Brain App validates and uses. Brain Web results must match direct model calls,
and generated data must remain outside Git.
