# Parameter rationale and reliability

The current model is still in a normalized-unit regime. That means many values
are effective rate constants rather than direct physical constants. To make this
honest and useful, the iterated model now keeps nominal values and uncertainty
ranges in `config/parameter_provenance.yaml`.

## Literature-informed anchors

- AAV is a widely used in vivo gene delivery vector, with serotype-dependent
  tropism and limited packaging capacity. AAV9-like delivery is often used as
  the conceptual basis for BBB-crossing CNS delivery.
  Sources:
  - https://en.wikipedia.org/wiki/Adeno-associated_virus
  - https://en.wikipedia.org/wiki/Gene_therapy

- BBB delivery is a major bottleneck; large biological cargos often require
  specialized crossing mechanisms, and transcytosis is one conceptual route.
  Sources:
  - https://en.wikipedia.org/wiki/Blood%E2%80%93brain_barrier
  - https://en.wikipedia.org/wiki/Transcytosis

- Mammalian mRNA and protein half-lives vary widely. Therefore, expression and
  degradation rates are represented as broad priors rather than fixed truths.
  Source:
  - https://en.wikipedia.org/wiki/Messenger_RNA

- APOBEC-family enzymes catalyze cytidine deamination and RNA/DNA editing, but
  engineered PUF-APOBEC activity is construct- and substrate-dependent. The
  editing parameters must therefore be fitted with amplicon sequencing.
  Sources:
  - https://en.wikipedia.org/wiki/APOBEC
  - https://en.wikipedia.org/wiki/APOBEC1
  - https://en.wikipedia.org/wiki/APOBEC3A

## Why broad priors are used

The exact parameters for our system depend on:

- AAV capsid or delivery vehicle actually used by the wet lab.
- Dose and route.
- Cell type.
- Promoter strength.
- PUF domain sequence and target RNA accessibility.
- APOBEC fusion design.
- RNA-seq or amplicon-seq off-target landscape.

Because these are not fixed yet, the best current practice is to:

1. Use biologically plausible ranges.
2. Run uncertainty analysis.
3. Identify which parameters dominate decisions.
4. Ask wet lab to measure those parameters first.
5. Refit the model after each experimental round.

This is exactly the DBTL loop implemented in the iterated model.
