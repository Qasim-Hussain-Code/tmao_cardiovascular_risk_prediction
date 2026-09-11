# Data dictionary: ST001420

Schema for the reanalysis cohort. The original study schema, which governs the
simulated cohort and any future dataset with clinical covariates, is in
[data_dictionary.md](data_dictionary.md). The two are not interchangeable and
neither loader accepts the other's file.

Source: Metabolomics Workbench study ST001420, project PR000974, deposited by
Cui and colleagues at the University of California, San Diego. Released 30
July 2020 under CC BY 4.0. Retrieved by `fetch_mwtab` in
`src/tmao_cvd/st001420.py` and written to `data/raw/`, which is not tracked.

## Measurement scale

**All metabolite values are relative peak areas in arbitrary instrument units.
They are not concentrations.**

No value in this dataset may be compared against a concentration threshold,
cut point, reference range or tertile boundary taken from any other study, and
no absolute level here may be called high or low relative to a published
cohort. Only comparisons internal to this dataset are admissible.

This is enforced rather than documented. `src/tmao_cvd/scales.py` declares the
scale as a typed value, `apply_concentration_threshold` raises
`MeasurementScaleError` when handed peak area data, and the loader verifies the
declared scale against the values on every load so that substituting a
calibrated file fails immediately. The reasoning is in section 4 of
[reanalysis_plan.md](reanalysis_plan.md).

## Columns

| Column | Type | Units | Description |
| --- | --- | --- | --- |
| `sample_id` (index) | string | | Deposit identifier, S1 to S750 |
| `recurrent_angina_9mo` | binary | | 1 if angina recurred within nine months of PCI, else 0. 210 events, 540 non-events |
| `sample_index` | integer | | Numeric part of `sample_id`. **Not a predictor.** Retained only for the acquisition order diagnostic in section 6 of the reanalysis plan, and excluded from `metabolite_columns` |
| 600 metabolite columns | numeric | relative peak area | Targeted LC-MS/MS, SCIEX QTRAP 6500+ |

### Pathway metabolites

These four form the models for question 1. All four are complete, with no
missing values in any of the 750 participants.

| Column | Role in the analysis |
| --- | --- |
| `Trimethylamine N-oxide` | The product under test |
| `Choline` | Dietary precursor, baseline model |
| `Betaine` | Choline oxidation product, baseline model |
| `Carnitine` | Dietary precursor, baseline model |

### Derived variables

Built by `add_log_pathway_features`. The log transform is monotone, so a
logged peak area is still an arbitrary unit and still carries no external
comparability.

| Column | Definition |
| --- | --- |
| `log_tmao` | natural log of `Trimethylamine N-oxide` |
| `log_choline` | natural log of `Choline` |
| `log_betaine` | natural log of `Betaine` |
| `log_carnitine` | natural log of `Carnitine` |

## Known structural problems

**Acquisition order is confounded with outcome.** Samples S1 to S210 are
exactly the 210 cases and S211 to S750 exactly the 540 controls. No batch
identifier, run order, injection sequence or acquisition date was deposited,
so instrument drift cannot be distinguished from biological signal or adjusted
for. Section 5 of the reanalysis plan sets out what this does and does not
permit a result from these data to mean. A test asserts the confounding so the
limitation cannot be forgotten.

**Panel missingness is high.** Across all 600 metabolites roughly 32.5 per
cent of cells are empty, and 195 metabolites are missing in more than a fifth
of participants. The four pathway metabolites are unaffected, so question 1 is
untouched by this. Question 2, which uses the whole panel, is not, and its
handling of missingness is specified in the plan rather than chosen after
seeing results.

**The endpoint is symptom driven.** Recurrent angina is reported and
adjudicated on symptoms, not a hard endpoint like myocardial infarction or
cardiovascular death. It is more susceptible to ascertainment differences than
the outcome the original study plan was written around.
