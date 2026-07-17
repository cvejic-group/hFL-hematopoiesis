# SCENIC+ analysis and Cell type-resolved eRegulon (CT-eRegulon) construction

- `scenicplus`: This folder contains code for running **SCENIC+** pipeline with **pycisTopic_polars_1xx** version
- `celltype_resolved_grn`: This folder contains code for converting the **raw SCENIC+ eRegulon set** into
**cell type-resolved eRegulons (CT-eRegulons)**. It corresponds to the *"Cell type-resolved eRegulon construction"* and
*"Cross-cell-type differential target analysis of CT-eRegulons"* sections in Methods.
- `dev_diff_eRegulon`: This folder contains code for developmental differential analysis of HSCs using 
linear mixed models (LMM) on +/+ activator eRegulons.
- `pycisTopic_distal_only`: This folder contains code for running pycisTopic on distal regions only.
- `Scripts`: 
  * `01.metaplot_Validation_eRegulon.R`: This script contains code for validating eRegulon regions using [public ChIp-seq data](https://ashpublications.org/blood/article/142/17/1448/497482/Genome-wide-transcription-factor-binding-maps).
  * `02.ERG_eRegulon_shRNA.R`: This script contains code for validating ERG eRegulon, specifically ERG-NKAIN2 link, using [ERG shRNA knocked down data](https://ashpublications.org/blood/article/138/16/1441/476046/Disruption-of-a-GATA2-TAL1-ERG-regulatory-circuit).

---

## Detailed explaination for how `celltype_resolved_grn` folder relates to the SCENIC+ pipeline

The CT-eRegulon construction is **not** part of the SCENIC+ workflow itself.
SCENIC+ provides the *input*; the scripts here implement a custom two-stage
filtering and pruning strategy on top of it.

```
                   ┌──────────────────────────────────────────────┐
   upstream        │  SCENIC+ (v1.0a1) + pycisTopic + pycistarget │
                   │  + GRNBoost2  →  317 high-confidence         │
                   │  eRegulons (TF–region–gene triplets, +/+,    │
                   │  +/-, -/+, -/-)                              │
                   └──────────────┬───────────────────────────────┘
                                  │
                                  │  raw outputs handed off as
                                  │  - eRegulon metadata table
                                  │ 
                                  ▼
                   ┌──────────────────────────────────────────────┐
   THIS FOLDER     │  Step 1 — one-vs-rest XGBoost + SHAP + RMT   │
   (custom code)   │           null  →  initial CT-eRegulons      │
                   │                                              │
                   │  Step 2 — pruning by cell-type DARs          │
                   │           →  final CT-eRegulons              │
                   │                                              │
                   │  Step 3 — Cross-cell-type analysis           |
                   │           (Jaccard + Walktrap modules)       │
                   └──────────────────────────────────────────────┘
```

### What is used as input from previous results

- **eRegulon metadata table** — the TF–region–gene triplet table that SCENIC+
  emits, which indicates the eRegulon/eGRN structure.
- **Cell-type peaks / DARs** called by ArchR `addReproduciblePeakSet` and
  `getMarkerFeatures`. These come from our annotation and preprocessing steps.

The scripts here do not modify any of those upstream objects.

### What is implemented from scratch in this folder

1. **One-vs-rest XGBoost classifier** on TF expression.
2. **SHAP-based TF ranking**, including the Gaussian/RMT null surrogate for
   significance.
3. **Tri-threshold filtering** (`logFC > 0`, `R_f^C > 2`,
   `ΔR_f > ΔR_cutoff`) to select representative TFs per cell type → filtered
   CT-eRegulons.
4. **DAR-based pruning**: drop target regions of an eRegulon that are
   not differentially accessible in the focal cell type, drop genes linked
   only through those regions, and discard CT-eRegulons that lose >90% of
   targets and are left with <5 genes → pruned CT-eRegulons.
5. **Cross-cell-type Jaccard analysis** of the pruned CT-eRegulons, plus
   **Walktrap module detection** on the resulting Jaccard similarity network.

---

## Reproducibility notes

- The XGBoost/SHAP step is stochastic (subsampling + tree training). Set
  the `seed` argument in the function prarmeters to
  reproduce the exact CT-eRegulon set reported in the manuscript.
- All downstream steps (DAR pruning, Jaccard, Walktrap) are deterministic
  given fixed inputs.
