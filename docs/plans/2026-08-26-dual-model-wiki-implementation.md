# Dual-Model iGEM Wiki Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a production-ready iGEM model Wiki, a Brain Delivery Web App, and a reproducible PUF-to-Brain candidate-panel workflow.

**Architecture:** A static React/Vite Wiki links to two independent FastAPI applications. PUF Atlas exports a stable CSV contract; Brain Delivery validates and aggregates that contract before running its existing v4 simulator.

**Tech Stack:** React, TypeScript, Vite, Vitest, Playwright, FastAPI, Jinja2, Plotly, pytest.

---

### Task 1: Candidate-panel contract

- Write failing PUF tests for unique-locus Top 100 export, conservative missing-value handling, metadata warnings, and Web downloads.
- Implement the exporter and add it to every completed scan.
- Run the focused tests, then all PUF tests.

### Task 2: Brain Delivery service layer

- Write failing tests for input validation, design application, panel validation/aggregation, single simulation, and dose scan.
- Implement typed service functions that call the existing v4 model without changing its equations.
- Verify Web/CLI numerical agreement.

### Task 3: Brain Delivery Web App

- Write failing route tests for home, simulation, optimization, uploads, downloads, and invalid input.
- Implement FastAPI routes, Jinja templates, local Plotly assets, and the scientific field-atlas styling.
- Add a documented `run_brain_app.py` entrypoint.

### Task 4: iGEM Wiki

- Scaffold the React/Vite site with tests for seven routes, configurable URLs, navigation, and required disclaimers.
- Implement the dual-track landing experience, shared scientific components, and English content sourced from repository documentation.
- Add self-hosted assets, responsive behavior, accessibility, and iGEM-compatible base paths.

### Task 5: Integration and verification

- Document local startup and environment variables; extend ignore rules for generated artifacts.
- Run pytest for both models, Python lint/type checks, Wiki lint/type/test/build, and asset/link audits.
- Run Playwright across desktop/mobile and the synthetic PUF-panel-Brain workflow.
