# Future Iterations

The v4 model is ready for wiki presentation, but it is also designed to grow.
This section can be used as the "Future Work" or "Next DBTL Cycle" part of the
wiki.

## Iteration A: HEK293T Digital Twin

Current full-chain model still includes AAV and BBB delivery. For immediate wet
lab, the highest-value next step is a HEK293T transfection-only mode:

```text
plasmid dose -> editor mRNA -> editor protein -> APOE editing
```

This would make the 24/48/72 h design screen directly match the wet-lab system.

## Iteration B: Real APOE Sequence Scanner

Current v4 uses a candidate design table. Next, the model should automatically:

- load APOE RNA sequence;
- enumerate candidate PUF binding sites;
- calculate editable C distance from PUF site;
- classify UC/CC/AC/GC context;
- predict RNA accessibility;
- search transcriptome-like sequences for off-target candidates.

This would replace hand-curated `puf_score` and `accessibility` values.

## Iteration C: Public REWIRE Benchmark

Use Han et al. 2022 public data to fit sequence-to-kinetics coefficients:

- target sequence;
- PUF repeat number;
- editing position;
- on-target editing percentage;
- bystander editing;
- transcriptome off-target.

This would turn current literature-informed rules into fitted benchmark rules.

## Iteration D: Wet-Lab Parameter Fitting

When our team has qPCR, Western blot, and amplicon sequencing, update:

- `k_tx`, `k_tl`, `k_deg_m`, `k_deg_p`;
- `k_cat_112`, `k_cat_158`;
- `local_bystander_per_112/158`;
- `k_on_puf_off`, `k_cat_puf_off`;
- `k_deaminase_bg`.

## Iteration E: Identifiability and Model Comparison

Before claiming exact parameters, add:

- parameter recovery;
- profile likelihood;
- prediction intervals;
- comparison between simple editing model and explicit binding/catalysis model.

This makes the model more rigorous and prevents overclaiming.

## Iteration F: In Vivo Translational Extension

If the project moves beyond HEK293T:

- calibrate AAV biodistribution with blood/liver/brain qPCR;
- separate brain vascular and parenchymal compartments;
- add liver/peripheral expression burden;
- test whether the chosen editor remains safe at systemic scale.

## Final Future Message

v4 is not the endpoint. It is a scaffold for repeated DBTL cycles: every new
experiment should either update a parameter, reduce uncertainty, or eliminate a
poor design.
