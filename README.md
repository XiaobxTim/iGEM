# REWIRE modeling workspace

This repository contains two independent scientific models and one static iGEM
Wiki that explains how they work together.

| Component | Directory | Default URL |
|---|---|---|
| PUF-OffTarget Atlas | `puf-offtarget-atlas/` | `http://127.0.0.1:8000` |
| Brain Delivery Model v4 | `brain_delivery_model/XiaobxTim-iGEM-REWIRE-iterated-v4/` | `http://127.0.0.1:8001` |
| Dual-model Wiki | `wiki/` | `http://127.0.0.1:5173` |

## Start all three interfaces

Use three terminals.

```bash
cd puf-offtarget-atlas
python -m pip install -e ".[dev]"
pufscan web
```

```bash
cd brain_delivery_model/XiaobxTim-iGEM-REWIRE-iterated-v4
python -m pip install -r requirements.txt
python run_brain_app.py
```

```bash
cd wiki
npm install
npm run dev
```

## PUF → Brain workflow

1. Complete a PUF-OffTarget Atlas scan.
2. Download `brain_candidate_panel.csv` and
   `brain_candidate_panel.metadata.json` from the result page.
3. Open the Brain Delivery App and attach the CSV in the optional candidate
   panel field.
4. Run a single simulation or 12-dose optimization, then download JSON or CSV.

The CSV transfer is manual and inspectable; the applications do not send data
to one another. Missing expression or accessibility scores are conservatively
set to `1.0` by the exporter and recorded in metadata.

Generated transcriptomes, scan results, model outputs, Wiki builds and package
dependencies are ignored by Git. Model outputs are literature-informed design
hypotheses, not clinically calibrated predictions.

## Integration verification

The automated cross-model test performs a synthetic PUF scan, parses its real
CSV export and runs the Brain model:

```bash
PYTHONDONTWRITEBYTECODE=1 conda run -n xbx_env pytest -q tests/test_cross_model_workflow.py
```

The browser smoke test expects the Wiki, Brain App and PUF App on ports 5173,
8001 and 8000 respectively:

```bash
conda run -n xbx_env python tests/browser_smoke.py
```
