# Plan: treatment-adjusted sensitivity analysis

Written and approved on 9 November 2026, before any sensitivity model was fitted. The author accepted
the proposed settings for decisions A to D ("go with your suggestions").
Reported as a post hoc analysis: it was added after the main results were known.

## Question
Do the Hallmark-programme associations with two-year PFS status, and the patient-level score and Cox
results, change materially when treatment is added to the adjustment set (UICC stage and HPV16 DNA/RNA
status)? The main results stay primary and unchanged.

## Cohorts (identical to the main analyses)
222 patients for the adjusted limma, ranked enrichment and sample-level score models; 250 patients with
120 PFS events for the Cox models. Treatment is recorded for all 270 patients (single-modality 78,
multi-modality 189, palliative 3); no imputation.

## Decision A: coding of treatment
Primary: multi-modality (reference) versus less than multi-modality (single-modality and palliative
combined). Fixed now; not changed after seeing results. Pre-specified check: repeat with the palliative
patients excluded (1 patient in the 222 cohort, 3 in the 250 cohort).

## Models (settings copied from the main analysis, plus one treatment term)
(a) limma: group + UICC stage + HPV16 DNA/RNA + treatment (222 patients); moderated t ranked; fgsea
    fgseaMultilevel, set.seed(123), minSize 15, maxSize 500, eps 0, same 50 Hallmark sets.
(b) Sample-level scores: linear model score ~ UICC stage + HPV16 DNA/RNA + treatment + group, 50 sets,
    Benjamini-Hochberg within the family.
(c) Cox: Surv(PFS) ~ standardised score + UICC stage + HPV16 DNA/RNA + treatment, Efron ties,
    cox.zph (km transform), Benjamini-Hochberg within the family of 50 (250 patients).

## Outcomes (fixed in advance)
P1. Of the 23 sets significant in the main ranked analysis, how many keep the sign of the NES.
P2. How many of those 23 remain at q < 0.05, and whether any other set becomes significant.
P3. Median and range of NES(sensitivity) / NES(main) across the 23 sets.
P4. Spearman correlation of gene-level t between main and sensitivity models; number of genes with q < 0.05.
P5. Sample-level and Cox models: for the nominal main models (9 score, 4 Cox) whether direction and
    p < 0.05 are kept; percentage change in the estimate (beta; log hazard ratio) for E2F Targets, MYC
    Targets V1, G2M Checkpoint and Estrogen Response Late; whether any model reaches q < 0.05.

## Decision B: how the result is described (descriptive benchmarks, not tests)
Robust: all 23 keep sign, at least 19 of 23 remain at q < 0.05, median NES ratio between 0.8 and 1.2.
Not robust: more than 2 sets reverse sign, or fewer than half (fewer than 12) remain at q < 0.05.
Partly attenuated: anything in between. The full results are reported whatever the outcome.

## Decision C: reporting
New S6 Table (all 50 sets from each model) plus one short Methods paragraph and one Results paragraph;
revised limitation sentence and abstract sentence.

## Not done
No other covariates (age, sex, smoking, alcohol, site), no alternative treatment codings after the
results are seen, no change to the main results, no claim of independence from treatment.

## Limits
Treatment is decided after the biopsy and reflects stage, site, fitness and biology (risk of
over-adjustment); the variable is crude; performance status and comorbidity are unavailable; treatment
is strongly associated with stage. The analysis shows whether findings move when treatment is added; it
cannot show independence from treatment.
