# Benchmark research: what strong iGEM modeling projects do well

This file summarizes design patterns borrowed from strong iGEM projects and
from our team lead's REWIRE RNA-editing pipeline structure.

## Reference projects and reusable strengths

### Heidelberg 2024: model as a digital twin

The Heidelberg 2024 model page presents DaVinci as an end-to-end digital twin:
structure prediction, molecular dynamics, long-range simulation, custom scoring
and wet-lab feedback. Their strongest modeling pattern is not any single tool,
but the closed loop: in silico design, wet-lab measurement, then model update.

Source: https://2024.igem.wiki/heidelberg/model

Useful practices imported into our model:

- Treat the model as an engineering decision tool, not just a plot generator.
- Make assumptions and parameter choices explicit.
- Use experimental readouts to calibrate or reject model assumptions.
- Report what the model cannot yet explain.

### Heidelberg 2024 Engineering Success: visible DBTL cycles

Their engineering page documents repeated Design-Build-Test-Learn cycles across
construct design, assays, burden reduction and model refinement. It shows
failed assays and changed designs, which makes the project credible.

Source: https://2024.igem.wiki/heidelberg/engineering

Useful practices imported into our model:

- Add DBTL records for each model iteration.
- Track failure modes instead of hiding them.
- Use wet-lab bottlenecks to decide which model feature comes next.

### Marburg 2021: model-driven biological design

Marburg 2021 used computational models to search sequence space, reduce design
complexity and propose synthetic regulatory parts. Their modeling page is strong
because it links algorithmic outputs to concrete biological design choices.

Source: https://2021.igem.org/Team:Marburg/Model

Useful practices imported into our model:

- Make model outputs actionable: recommended assays, parameter priorities,
  dose windows and safety trade-offs.
- Keep a clear interface between dry lab and wet lab.

### REWIRE RNA-editing pipeline: reproducible dry-lab repository structure

The reference repository is especially strong structurally. It separates
pipeline code, DBTL records, decision logs, failure logs, wiki material, result
summaries and tests.

Source: https://github.com/pdx12320/REWIRE-RNA-editing-pipeline

Useful practices imported into our model:

- Add `dbtl/`, `docs/`, `wetlab/`, `reports/` and script-level entry points.
- Freeze output summaries in CSV/Markdown.
- Provide a wet-lab template so future data can enter the model cleanly.

## What this means for our AAV-PUF-APOBEC model

The old combined model had a useful mechanistic chain:

```text
AAV delivery -> brain expression -> PUF-APOBEC RNA editing
```

The iterated model adds the missing iGEM-quality layer:

```text
assumption register
parameter provenance
wet-lab observation interface
uncertainty scan
model-informed next-experiment recommendation
DBTL records
```

This makes the model easier to defend on a wiki: judges can see not only the
ODE equations, but also why parameters were chosen, which wet-lab data would
test them, and what the next iteration should measure.
