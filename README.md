# TMAO, recurrent angina, and what a public deposit can support

[![tests](https://github.com/Qasim-Hussain-Code/tmao_cardiovascular_risk_prediction/actions/workflows/tests.yml/badge.svg)](https://github.com/Qasim-Hussain-Code/tmao_cardiovascular_risk_prediction/actions/workflows/tests.yml)

## 1. What this chapter set out to answer?

Trimethylamine N-oxide is produced when gut bacteria convert dietary choline,
phosphatidylcholine and L-carnitine into trimethylamine, which the liver then
oxidises through flavin-containing monooxygenase 3. Higher circulating
concentrations have been associated with incident cardiovascular events in
several prospective cohorts.

Association is not predictive usefulness. A marker earns a place in a risk
model only if it improves discrimination, calibration or clinical decision
making beyond the variables a clinician already holds. The original question,
fixed in [docs/analysis_plan.md](docs/analysis_plan.md), was therefore narrow:
does plasma TMAO add predictive value to a conventional cardiovascular risk
model for incident major adverse cardiovascular events in primary prevention?

TMAO is a pointed case because it is cleared renally. Impaired renal function
raises circulating TMAO and independently raises cardiovascular risk, so any
honest assessment must hold renal function fixed before asking what TMAO
contributes.

## 2. Why that question could not be answered, and what the search found

No public dataset supports it. The cohorts that could, chiefly the
Multi-Ethnic Study of Atherosclerosis and the Cardiovascular Health Study,
hold their TMAO measurements under controlled access requiring institutional
affiliation, ethical approval and a data use agreement. TMAO is not part of
the standard BioLINCC release for those cohorts, so even an approved request
for the parent study would not return the one variable this chapter is about.

The first search was unsuccessful for a reason worth recording, because it was
a methodological error rather than bad luck. I searched repository *study
titles* for trimethylamine and TMAO, which returned three studies, only one of
them human and that one with 36 participants and no outcome. Searching instead
by *measured metabolite* returned 98 human blood studies. The distinction is
that a metabolite appears in the panel of studies whose titles never mention
it, which is the ordinary case for an untargeted or broad targeted assay.

That second search located a public Metabolomics Workbench deposit, accession
ST001420, released under CC BY 4.0. The depositing group is not named anywhere
in this repository, and nothing here is a comment on them.

## 3. What the deposit can and cannot answer

| Property | Value |
| --- | --- |
| Design | Prospective observational cohort, sampled after percutaneous coronary intervention |
| Participants | 750 |
| Events | 210 recurrent angina by nine months |
| Panel | 600 metabolites by targeted mass spectrometry |
| Licence | CC BY 4.0 |

**It cannot answer the original question.** The deposit carries no
participant-level clinical covariates. Age, sex, smoking and lipids exist only
as aggregate summaries in the study description, so there is no conventional
risk model for TMAO to be incremental to.

**It can answer a narrower one**, fixed in advance in
[docs/reanalysis_plan.md](docs/reanalysis_plan.md) and committed before the
code it governs: does TMAO carry information beyond its own dietary
precursors, in a secondary prevention cohort of post-PCI patients with stable
angina, of Asian ethnicity, uniformly on dual antiplatelet therapy, predicting
recurrent angina within nine months?

That is a different claim, not a weaker version of the original. Nothing
learned here transfers to primary prevention, to hard cardiovascular
endpoints, or to other populations.

Two constraints bound every result. Values are relative peak areas in
arbitrary instrument units rather than concentrations, so no published cut
point may be applied to them; this is enforced by an exception in
[src/tmao_cvd/scales.py](src/tmao_cvd/scales.py), not by a comment. And sample
identifiers are perfectly ordered by outcome, with no batch identifier, run
order or acquisition date deposited.

### The analysis at a glance

![Workflow of the reanalysis: the public deposit and what it omits on the left,
the six pipeline stages through validation, the three pre-specified questions
and reporting in the centre, and the summary file, verdict report, per question
tables and audit result on the right. The post hoc diagnostics branch off
question two.](docs/figures/workflow_diagram.jpeg)

The post hoc diagnostics hang off question two as a branch rather than sitting
in the numbered sequence. That is deliberate, and it is the honest shape of
what happened: they were triggered by a result, not planned before one existed.
The distinction is preserved in the code layout and in the output tables as
well as in this diagram.

## 4. The pre-specified questions and their verdicts

Two substantive questions and one diagnostic, each with its verdict decided by
thresholds fixed before any estimate existed.

### Question 1, primary. Verdict: POSITIVE

Baseline model is log choline, log betaine and log carnitine. The extended
model adds log TMAO and nothing else.

| Model | Area under the curve | Interval |
| --- | --- | --- |
| Baseline, precursors only | 0.8475 | 0.8132 to 0.8818 |
| Extended, plus log TMAO | 0.8731 | 0.8402 to 0.906 |
| Difference | 0.0256 | 0.0075 to 0.0436 |

The interval on the difference excludes zero, with a paired DeLong p of
0.00547. Integrated discrimination improvement is 0.0669, interval 0.0454 to
0.0882. Category-free net reclassification improvement totals 0.6571, interval
0.4973 to 0.8016. The extended model's decision curve lies above the
baseline's across a proportion 0.9821 of the pre-specified threshold range.

### Question 2, replication. Verdict: PARTIAL

An L2 penalised model over the whole panel, with the penalty chosen inside
each training fold, evaluated strictly out of fold.

| Metric | Out of fold |
| --- | --- |
| Accuracy | 0.9973 |
| Sensitivity | 1.0 |
| Specificity | 0.9963 |
| Area under the curve | 1.0 |
| Calibration slope | 2.9711 |
| Calibration intercept | 0.0218 |

A figure of above 0.890 on accuracy, sensitivity and specificity has
previously been reported in connection with this deposit. Discrimination here
sits above the pre-specified reproduction floor, but the calibration slope of
2.9711 falls far outside the pre-specified band, so the reproduction criterion
is not met and the verdict is PARTIAL. An area under the curve of 1.0 on 750
participants is not a performance result to be reported and moved past. It is
a signal that something structural separates the two groups, which is what
question 3 and the post hoc diagnostics then examined.

### Question 3, diagnostic. Verdict: ORDER STRUCTURE DETECTED

Spearman correlation between sample index and log peak area, computed within
each outcome block separately so the test is blind to the between-group
difference. Under the null of no order structure the p values are uniform, so
the expected proportion below 0.05 is 0.05.

| Block | Participants | Metabolites tested | Proportion at p below 0.05 | Kolmogorov-Smirnov p |
| --- | --- | --- | --- | --- |
| Cases | 210 | 405 | 0.1679 | 5.88e-08 |
| Controls | 540 | 405 | 0.1556 | 8.3e-11 |

Both blocks breach both pre-specified limits.

### Post hoc diagnostics

Not pre-specified. Prompted by the question 2 result, and labelled as post hoc
wherever they appear.

| Diagnostic | Value | Under the null |
| --- | --- | --- |
| Median univariate area under the curve across the panel | 0.6852 | 0.5 |
| Proportion of metabolites above 0.6 | 0.7654 | |
| Proportion above 0.8 | 0.2444 | |
| Proportion above 0.9 | 0.0716 | |
| Total log intensity, direction-free area | 0.7289 | |

Missingness was the first hypothesis for the perfect separation and is not the
mechanism. Of the panel, 405 metabolites are fully observed and 195 are absent
for every participant, and the maximum difference in missingness between the
groups is 0.0.

## 5. The question 1 inversion

Two statements are both true, and reporting either without the other misleads.

**The statistic is positive exactly as pre-specified.** The interval on the
difference excludes zero, the decision curves separate, and the verdict fired
POSITIVE under rules fixed before any estimate existed. It has not been
revised after the fact.

**It is not evidence about TMAO biology.** TMAO has a univariate area under
the curve of 0.7404, which places it at the 64.2nd percentile of the 405 fully
observed metabolites, with 144 ranking above it. Among those ranking above it
is betaine, at 0.8427, which is one of its own precursors and therefore part
of the baseline model it was tested against. The panel median is 0.6852
against a null of 0.5, and a proportion 0.7654 of metabolites discriminate
above 0.6. An incremental contribution measured inside a panel displaced that
far is not evidence about the marker.

The plan originally argued the opposite, and that reasoning was wrong. Section
5 held that question 1 would be relatively robust to acquisition artefact,
because a drift affecting the whole panel similarly would largely cancel in a
comparison between metabolites from the same run. The premise fails: the panel
is affected heterogeneously, with univariate discrimination ranging from
chance to well above 0.9. Under heterogeneous drift one metabolite outperforms
another by being more sensitive to it, not by being more informative. The
correction is recorded in the departures section of the plan rather than
absorbed silently, and the verdict is reported unchanged with the
qualification attached.

This coupling is enforced in code. A single constructor builds the question 1
statement, the panel context is a required argument with no default, and a
test asserts no second code path can emit the difference without it.

## 6. The question 3 finding

The deposit contains no batch identifier, no run order, no injection sequence
and no acquisition date. Sample identifiers are perfectly ordered by outcome.
Consequently, if the samples were acquired in the order deposited, instrument
drift over the run is aligned with the outcome and cannot be distinguished
from biological signal; and because no batch metadata exists, that possibility
can neither be tested directly nor adjusted for.

This is a limit of what this particular dataset can tell us, not a flaw unique
to it. Like many public deposits, it does not include the laboratory processing
order alongside the results. Without that, there is no way to fully rule out
that something about how the two groups were handled, rather than their
biology, is doing some of the work. Batch metadata is routinely absent from
deposits of this vintage, and the finding here is about what a reader can
verify from the public record rather than about what was done in the
laboratory. It is not a claim about anyone who produced, analysed or deposited
the data.

The consequence is symmetrical, and worth stating plainly in both directions.
Order structure being detected does not establish that the signal is
artefactual. Had it not been detected, that would not have established the
signal was biological either, since the deposited ordering need not be the
acquisition ordering. The deposit does not contain what would settle it either
way. That is the finding.

## 7. What this chapter found

This chapter's finding is about how to check a claim, not about TMAO.

## The chapter as published

The seven posts that make up this chapter are in [posts/](posts/), committed
verbatim as published.

| Post | Subject |
| --- | --- |
| [day47.txt](posts/day47.txt) | The pathway, and the question the chapter opened with |
| [day48.txt](posts/day48.txt) | Why the gene does not follow the family tree |
| [day49.txt](posts/day49.txt) | The rule, written before anything was downloaded |
| [day50.txt](posts/day50.txt) | Four ways the plan could be wrong |
| [day51.txt](posts/day51.txt) | Why a support vector machine, and what the kernel does |
| [day52.txt](posts/day52.txt) | The pivot, stated publicly, and the narrower question |
| [day53.txt](posts/day53.txt) | The two results, and why the second is the one to trust |

Days 47 to 51 set out a study of gut metagenomic gene abundances that the
public data could not support. Day 52 states the pivot in public rather than
absorbing it quietly, and the work in this repository implements what that post
describes. The original question is not answered here and remains open.

## Reproducing

Python 3.10 or later, under either resolver.

```bash
conda env create -f environment.yml && conda activate tmao-cvd   # conda
python -m pip install -r requirements.txt                        # pip

python scripts/run_reanalysis.py    # downloads the deposit, runs all three questions
python scripts/audit_readme.py      # verifies every figure above against results/metrics
python -m pytest tests/             # test suite
```

The deposit is downloaded to `data/raw/`, which is not tracked. Generated
figures and tables are not tracked either, since they are reproducible from
the code and the recorded seed. `results/metrics/` is the exception and is
tracked, because it is the audit trail for this document.

## Layout

```
.
├── .github/workflows/      test suite run on every push
├── config/                 analysis settings kept outside the code
├── data/                   cohort files, never committed
├── docs/
│   ├── analysis_plan.md            the original pre-specified plan
│   ├── reanalysis_plan.md          the plan governing this chapter
│   ├── data_dictionary.md          schema for the original question
│   ├── data_dictionary_st001420.md schema for the deposit
│   └── figures/                    the workflow diagram
├── results/metrics/        machine readable summary, tracked
├── scripts/                entry points and the README audit
├── src/tmao_cvd/           the package
└── tests/                  unit, integration and guard tests
```

## Related reading

Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a
multivariable prediction model for individual prognosis or diagnosis (TRIPOD).
Annals of Internal Medicine. 2015;162(1):55-63.

DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more
correlated receiver operating characteristic curves: a nonparametric approach.
Biometrics. 1988;44(3):837-845.

Hlatky MA, Greenland P, Arnett DK, et al. Criteria for evaluation of novel
markers of cardiovascular risk. Circulation. 2009;119(17):2408-2416.

Pencina MJ, D'Agostino RB Sr, Steyerberg EW. Extensions of net
reclassification improvement calculations to measure usefulness of new
biomarkers. Statistics in Medicine. 2011;30(1):11-21.

Steyerberg EW, Vickers AJ, Cook NR, et al. Assessing the performance of
prediction models: a framework for traditional and novel measures.
Epidemiology. 2010;21(1):128-138.

Sun X, Xu W. Fast implementation of DeLong's algorithm for comparing the areas
under correlated receiver operating characteristic curves. IEEE Signal
Processing Letters. 2014;21(11):1389-1393.

Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating
prediction models. Medical Decision Making. 2006;26(6):565-574.

## Licence

MIT. See [LICENSE](LICENSE).
