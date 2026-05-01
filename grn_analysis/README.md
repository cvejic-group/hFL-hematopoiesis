# Cell type-resolved eRegulon (CT-eRegulon) construction

This folder contains the code that converts the **raw SCENIC+ eRegulon set** into
**cell type-resolved eRegulons (CT-eRegulons)** for the foetal-liver multiome
atlas. It corresponds to the *"Cell type-resolved eRegulon construction"* and
*"Cross-cell-type differential target analysis of CT-eRegulons"* sections of the
Methods, and produces the regulons summarised in **Fig. 1d–e** and
**Supplementary Fig. 5**.

---

## 1. How this folder relates to the SCENIC+ pipeline

The CT-eRegulon construction is **not** part of the SCENIC+ workflow itself.
SCENIC+ provides the *input*; the scripts here implement a custom two-stage
filtering and pruning strategy on top of it.

```
                   ┌──────────────────────────────────────────────┐
   upstream        │  SCENIC+ (v1.0a1) + pycisTopic + pycistarget │
   (not in this    │  + GRNBoost2  →  317 high-confidence         │
    folder)        │  eRegulons (TF–region–gene triplets, +/+,    │
                   │  +/-, -/+, -/-)  +  AUCell scores            │
                   └──────────────┬───────────────────────────────┘
                                  │
                                  │  raw outputs handed off as
                                  │  - eRegulon metadata table
                                  │  - gene-/region-based AUCell matrices
                                  ▼
                   ┌──────────────────────────────────────────────┐
   THIS FOLDER     │  Step 1 — one-vs-rest XGBoost + SHAP + RMT   │
   (custom code)   │           null  →  initial CT-eRegulons      │
                   │                                              │
                   │  Step 2 — pruning by cell-type DARs          │
                   │           →  final CT-eRegulons              │
                   │                                              │
                   │  Cross-cell-type analysis (Jaccard +         │
                   │  Walktrap modules)                           │
                   └──────────────────────────────────────────────┘
```

### What is reused from SCENIC+ (essentially copy/paste of upstream outputs)

These objects are produced upstream and consumed *as-is* by the scripts in this
folder:

- **eRegulon metadata table** — the TF–region–gene triplet table that SCENIC+
  emits (the same content as `eRegulon_Raw_Metadata.xlsx` in the project, with
  one sheet per regulation class: `++`, `+-`, `-+`, `--`). Only the `+/+`
  activator subset is used downstream (see Methods, Step 1).
- **Gene-based and region-based AUCell matrices** computed by SCENIC+ /
  pycisTopic on the 575,177-peak consensus set.
- **Cell-type peaks / DARs** called by ArchR `addReproduciblePeakSet` and
  `getMarkerFeatures`. These come from `../03.archr/` (see top-level repo
  `README`) and are used only for the Step 2 pruning.

The scripts here do not modify any of those upstream objects — they only read
them.

### What is implemented from scratch in this folder

Everything below the dashed line in the diagram. None of it is part of SCENIC+:

1. **One-vs-rest XGBoost classifier** on TF expression (Step 1).
2. **SHAP-based TF ranking**, including the Gaussian/RMT null surrogate for
   significance (`E_null` with i.i.d. N(0,1) entries) and the
   `R_f^C`, `R_f^{non-C}`, `ΔR_f` statistics.
3. **Tri-threshold filtering** (`logFC > 0`, `R_f^C > 2`,
   `ΔR_f > ΔR_cutoff`) to select representative TFs per cell type → initial
   CT-eRegulons.
4. **DAR-based pruning** (Step 2): drop target regions of an eRegulon that are
   not differentially accessible in the focal cell type, drop genes linked
   only through those regions, and discard CT-eRegulons that lose >90% of
   targets and are left with <5 genes.
5. **Cross-cell-type Jaccard analysis** of pruned target sets, plus
   **Walktrap module detection** on the resulting weighted similarity network.

If you are looking for SCENIC+ behaviour, look upstream (in this repo:
`grn_analysis/scenicplus/` — adjust path if different); if you are looking for
how a TF gets assigned to a specific cell type, the logic is here.

---

## 2. Inputs expected by these scripts

| Object                              | Origin                              | Used in       |
| ----------------------------------- | ----------------------------------- | ------------- |
| eRegulon metadata table             | SCENIC+ (`+/+` sheet)               | Step 1, Step 2 |
| TF expression matrix (cells × TFs)  | Seurat / Scanpy log-normalised RNA  | Step 1        |
| Cell-type labels                    | WNN annotation (`01.wnn_anno`)      | Step 1, Step 2 |
| Cell-type DAR sets                  | ArchR (`03.archr`)                  | Step 2        |
| Gene/region AUCell matrices         | SCENIC+ / pycisTopic                | sanity checks |

Replace the paths at the top of each script with your own locations.

---

## 3. Pipeline overview (file by file)

> Adjust the file names below to match your actual scripts.

- `step1_xgboost_shap.<py|ipynb>` — trains the per-cell-type one-vs-rest
  XGBoost classifier on TF expression, computes SHAP values, fits the
  Gaussian null, computes `R_f^C / R_f^{non-C} / ΔR_f` and `logFC`, applies
  the tri-threshold, and writes the initial CT-eRegulon table.
- `step2_dar_pruning.<py|ipynb>` — for each (cell type, TF) pair retained
  after Step 1, intersects target regions with that cell type's DARs, drops
  orphan target genes, and applies the size/retention filter
  (≥5 genes or <90% loss). Writes the **final CT-eRegulon table**.
- `crosscelltype_jaccard_walktrap.<py|ipynb>` — for every TF retained in ≥2
  cell types, builds the pairwise Jaccard matrix of pruned targets and runs
  Walktrap (`igraph`) to recover regulatory modules. Used for
  Supplementary Fig. 5b–d.
- `utils.py` — shared helpers (subsampling, RMT null generation, threshold
  selection).

---

## 4. Outputs

- `CT_eRegulons_initial.tsv` — TF × cell-type table of representative
  eRegulons after Step 1.
- `CT_eRegulons_final.tsv` — same after DAR pruning (Step 2). This is the
  table cited as "CT-eRegulons" in the main text and figures.
- `jaccard_matrices/` — one Jaccard matrix per recurrent TF.
- `walktrap_modules.tsv` — module assignment per (TF, cell type).

---

## 5. Key parameters (defaults match the manuscript)

- XGBoost: one-vs-rest binary classifier, balanced random subsampling per
  iteration (number of iterations: see script), default tree booster.
- SHAP: per-cell SHAP values averaged within / outside the focal cell type;
  ε = 1e-6 in `R_f` denominators.
- Null: Gaussian surrogate `E_null ~ N(0,1)` with the same shape and feature
  names as the TF expression matrix.
- Tri-threshold: `logFC > 0`, `R_f^C > 2`, `ΔR_f > ΔR_cutoff`
  (`ΔR_cutoff` chosen as described in the script — adjust here if you change
  it).
- Pruning retention: keep CT-eRegulon if pruned set has ≥5 target genes or
  ≤90% of original targets were removed.
- Walktrap: run on the Jaccard-weighted network per TF; default step length.

---

## 6. Citation

If you use the CT-eRegulon code, please cite this manuscript and the SCENIC+
suite (the upstream eGRN inference):

- Bravo González-Blas *et al.*, *SCENIC+: single-cell multiomic inference of
  enhancers and gene regulatory networks*, Nat. Methods (2023).
- pycisTopic, pycistarget, GRNBoost2 — see SCENIC+ documentation.

---

## 7. Reproducibility notes

- The XGBoost/SHAP step is stochastic (subsampling + tree training). Set
  the `random_state` / `seed` arguments at the top of the Step 1 script to
  reproduce the exact CT-eRegulon set reported in the manuscript.
- All downstream steps (DAR pruning, Jaccard, Walktrap) are deterministic
  given fixed inputs.
