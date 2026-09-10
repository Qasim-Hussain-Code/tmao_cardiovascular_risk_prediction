# Data dictionary

Every dataset entering the pipeline, whether observed or simulated, must
carry these columns under these names and in these units. The schema is
enforced in code by `validate_cohort` in `src/tmao_cvd/data.py`, and the
ranges below are the ones that function checks.

The ranges are wide. They exist to catch unit errors and sentinel codes such
as 999, not to exclude values that are unusual but real. A participant whose
measurement falls outside a range should be investigated rather than quietly
dropped, which is why validation raises rather than filters.

| Column | Type | Units | Plausible range | Description |
| --- | --- | --- | --- | --- |
| `participant_id` | string | | | Unique identifier. Must not repeat. |
| `age_years` | numeric | years | 18 to 110 | Age at blood sampling. |
| `sex_male` | binary | | 0 or 1 | 1 for male, 0 for female, as recorded at enrolment. |
| `current_smoker` | binary | | 0 or 1 | Current smoking at baseline. |
| `diabetes` | binary | | 0 or 1 | Diagnosed diabetes mellitus at baseline. |
| `bmi_kg_m2` | numeric | kg per square metre | 10 to 70 | Body mass index. |
| `systolic_bp_mmhg` | numeric | mmHg | 60 to 260 | Seated systolic blood pressure. |
| `total_cholesterol_mmol_l` | numeric | mmol per litre | 1 to 15 | Total cholesterol. |
| `hdl_cholesterol_mmol_l` | numeric | mmol per litre | 0.2 to 5 | High density lipoprotein cholesterol. |
| `egfr_ml_min_1_73m2` | numeric | mL per min per 1.73 square metres | 5 to 200 | Estimated glomerular filtration rate. |
| `hs_crp_mg_l` | numeric | mg per litre | 0.01 to 200 | High sensitivity C reactive protein. |
| `tmao_umol_l` | numeric | micromoles per litre | 0.05 to 500 | Fasting plasma trimethylamine N-oxide. |
| `mace_3yr` | binary | | 0 or 1 | Major adverse cardiovascular event within three years of sampling. |

## Derived variables

Two variables are constructed by `add_derived_features` and are not expected
in the input file.

| Column | Definition | Reason |
| --- | --- | --- |
| `log_tmao` | natural log of `tmao_umol_l` | TMAO is strongly right skewed. On the untransformed scale a small number of very high values dominates the fit of a linear term. |
| `log_hs_crp` | natural log of `hs_crp_mg_l` | Same reasoning, and the log scale is the convention in the inflammatory marker literature. |

## Units worth checking before an export

Cholesterol is recorded here in mmol per litre. A file reporting mg per
decilitre will pass a careless eye and fail validation, which is the
intended behaviour. The conversion for total and HDL cholesterol is to
divide mg per decilitre by 38.67.

TMAO is recorded in micromoles per litre. Reported medians in adult cohorts
usually fall between about 2 and 5, so a column with a median in the
hundreds is most likely in nanomoles per litre.
