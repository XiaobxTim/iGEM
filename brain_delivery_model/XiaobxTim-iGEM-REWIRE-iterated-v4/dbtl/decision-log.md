# Decision log

| date | decision | reason |
|---|---|---|
| 2026-07-10 | Create a separate iterated version | Keeps old combined model intact and makes the new version easy to compare |
| 2026-07-10 | Use `P_brain` as free active editor | This is the cleanest interface between AAV expression and PUF-APOBEC editing |
| 2026-07-10 | Add wet-lab observation bridge before fitting | Fitting without a documented assay interface would be hard to defend |
| 2026-07-10 | Use broad parameter priors | Current wet-lab system is not calibrated, so exact-looking parameters would be misleading |
| 2026-07-10 | Recommend next experiments from uncertainty scan | Makes the model useful for DBTL planning |
| 2026-07-11 | Add REWIRE-style off-target candidate panel | Off-target risk should come from candidate sites rather than a magic aggregate number |
| 2026-07-11 | Add transparent random-search fitting | Early wet-lab data will be sparse, so an explainable fitter is better than a black-box optimizer |
