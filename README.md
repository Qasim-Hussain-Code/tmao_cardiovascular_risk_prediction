# TMAO and cardiovascular risk prediction

Analysis code for a single question: does fasting plasma trimethylamine
N-oxide (TMAO) add predictive value to a cardiovascular risk model that
already contains the variables a clinician has to hand?

## Background

TMAO is a small molecule of partly microbial origin. Gut bacteria convert
dietary choline, phosphatidylcholine and L-carnitine into trimethylamine,
which the liver then oxidises to TMAO through flavin containing
monooxygenase 3. Work from Hazen's group established the pathway and linked
it to atherosclerosis in animal models and to incident cardiovascular events
in clinical cohorts (Wang et al., 2011; Tang et al., 2013; Koeth et al.,
2013).

That an association exists is now reasonably well supported. Whether it is
useful for prediction is a separate question, and a harder one. A marker can
be strongly and reproducibly associated with an outcome and still add almost
nothing to a model containing age and kidney function. TMAO is a pointed
case, because it is cleared renally: impaired renal function raises
circulating TMAO and independently raises cardiovascular risk (Tang et al.,
2015). Any honest assessment has to hold renal function fixed before asking
what TMAO contributes.

This repository implements that assessment.

## Current status

There are no cohort data here, and none will be committed. What the code
does today is run end to end against a simulated cohort so that the analysis
can be written, tested and reviewed before a real dataset arrives.

This matters for how the output should be read. Running the pipeline as it
stands produces plausible looking tables and figures, and every one of them
is labelled `SIMULATED DATA, not an observed cohort` in the file itself. The
simulator plants a TMAO effect of its own choosing and the analysis recovers
it, which demonstrates that the code works. It demonstrates nothing whatever
about TMAO.

To run the analysis on real data, put the export at `data/raw/cohort.csv`
following `docs/data_dictionary.md`. The pipeline uses it automatically and
the provenance labels change accordingly.

## The comparison

Two nested models, differing by exactly one term.

- **Baseline.** Age, sex, current smoking, diabetes, body mass index,
  systolic blood pressure, total cholesterol, HDL cholesterol, estimated
  glomerular filtration rate, log C reactive protein.
- **Extended.** The baseline model plus log TMAO.

Keeping everything else identical is what makes the difference between them
readable as the contribution of that one term. The full specification, fixed
in advance, is in [docs/analysis_plan.md](docs/analysis_plan.md).

## How performance is judged

No single number settles whether a marker is worth measuring, so four
questions are asked separately, following the framework of Steyerberg and
colleagues (2010).

| Question | Measure |
| --- | --- |
| Does the model rank cases above controls? | Area under the ROC curve, compared between models by the paired DeLong test |
| Are the predicted risks numerically right? | Calibration intercept and slope, Brier score, Brier skill score |
| Do individuals move to more accurate risks? | Integrated discrimination improvement, category free net reclassification improvement, with bootstrap intervals |
| Would decisions actually improve? | Net benefit across risk thresholds, against treating everyone and treating nobody |

All performance figures come from out of fold predictions under repeated
stratified cross validation. Apparent performance is never reported, because
it is optimistic by an amount that grows with the number of candidate terms,
which is exactly the situation being studied.

The reclassification measures are included because the biomarker literature
expects them, not because they are trusted. Both are sensitive to
miscalibration and both reward predictions that are merely more extreme. The
decision curve carries more weight in interpretation.

## Layout

```
.
├── config/                 analysis settings kept outside the code
├── data/
│   ├── raw/                cohort export, never committed
│   ├── interim/            intermediate files, disposable
│   └── processed/          analysis ready tables, disposable
├── docs/
│   ├── analysis_plan.md    the prespecified statistical analysis plan
│   └── data_dictionary.md  column names, units and accepted ranges
├── results/
│   ├── figures/            generated figures
│   └── tables/             generated tables and the run manifest
├── scripts/
│   └── run_analysis.py     command line entry point
├── src/tmao_cvd/
│   ├── config.py           paths, seeds and resampling settings
│   ├── simulate.py         synthetic cohort generator
│   ├── data.py             loading and schema validation
│   ├── features.py         transformations and feature sets
│   ├── models.py           model specifications and resampling
│   ├── evaluate.py         discrimination, calibration, utility
│   ├── figures.py          the three reporting figures
│   └── pipeline.py         orchestration
└── tests/                  unit and end to end tests
```

## Installation

Python 3.10 or later.

```bash
python -m pip install -r requirements.txt
```

## Running

```bash
python scripts/run_analysis.py          # full analysis
python scripts/run_analysis.py --help   # seed, repeats, bootstrap, cohort path
python -m pytest tests/                 # test suite
```

A run takes about fifteen seconds on a laptop and writes four tables, three
figures and a manifest recording the seed, the software versions and the
data provenance behind them.

## What a simulated run should reproduce

With the default seed, the pipeline reports a baseline area under the curve
of 0.753 and an extended area of 0.762, a difference of 0.009 with a DeLong
p value of 0.013. These figures are useful only for checking that an
installation behaves as expected. They are properties of the simulator.

## Testing

Thirty tests cover the statistics and the joins between stages. The two
worth singling out:

The area under the curve and its variance are implemented here rather than
taken from a library, because the paired DeLong comparison is not available
in scikit-learn. The implementation is checked against `roc_auc_score` for
agreement to nine decimal places, and separately on the tied case, since
midranks only matter when predictions tie.

`tests/test_recovery.py` fits the whole chain to a large synthetic cohort
and asserts that every planted log odds ratio comes back within 0.05. This
catches the class of error that unit tests miss: a transformation applied
twice, a feature set assembled from the wrong columns, a scaler fitted on
the wrong axis. Each of those leaves the individual functions correct and
the analysis wrong.

## Reproducibility

Every stochastic step draws from a single seed recorded in
`src/tmao_cvd/config.py` and written into `results/tables/run_manifest.json`
alongside the library versions and the data provenance. Generated outputs
are not committed, since they are reproducible from the code and the seed,
and a stale figure in version control is worse than no figure at all.

## Limitations

The outcome is modelled as a binary indicator at three years. If a real
cohort carries meaningful loss to follow up this is the wrong tool, and the
analysis should move to a time to event model with administrative censoring.

Missing values are imputed at the median within each training fold, which is
adequate only for a small proportion of missingness unrelated to the
outcome. Beyond roughly ten per cent, median imputation understates the
uncertainty and multiple imputation is needed.

The DeLong confidence interval uses the normal approximation on the area
scale, which becomes conservative as the area approaches one.

Cross validation estimates how well this modelling procedure performs on
data drawn from the same population. It says nothing about transportability
to a different setting, which requires external validation.

## References

Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a
multivariable prediction model for individual prognosis or diagnosis
(TRIPOD). Annals of Internal Medicine. 2015;162(1):55-63.

DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or
more correlated receiver operating characteristic curves: a nonparametric
approach. Biometrics. 1988;44(3):837-845.

Hlatky MA, Greenland P, Arnett DK, et al. Criteria for evaluation of novel
markers of cardiovascular risk: a scientific statement from the American
Heart Association. Circulation. 2009;119(17):2408-2416.

Koeth RA, Wang Z, Levison BS, et al. Intestinal microbiota metabolism of
L-carnitine, a nutrient in red meat, promotes atherosclerosis. Nature
Medicine. 2013;19(5):576-585.

Pencina MJ, D'Agostino RB Sr, D'Agostino RB Jr, Vasan RS. Evaluating the
added predictive ability of a new marker: from area under the ROC curve to
reclassification and beyond. Statistics in Medicine. 2008;27(2):157-172.

Pencina MJ, D'Agostino RB Sr, Steyerberg EW. Extensions of net
reclassification improvement calculations to measure usefulness of new
biomarkers. Statistics in Medicine. 2011;30(1):11-21.

Steyerberg EW, Vickers AJ, Cook NR, et al. Assessing the performance of
prediction models: a framework for traditional and novel measures.
Epidemiology. 2010;21(1):128-138.

Sun X, Xu W. Fast implementation of DeLong's algorithm for comparing the
areas under correlated receiver operating characteristic curves. IEEE Signal
Processing Letters. 2014;21(11):1389-1393.

Tang WHW, Wang Z, Levison BS, et al. Intestinal microbial metabolism of
phosphatidylcholine and cardiovascular risk. New England Journal of
Medicine. 2013;368(17):1575-1584.

Tang WHW, Wang Z, Kennedy DJ, et al. Gut microbiota-dependent trimethylamine
N-oxide (TMAO) pathway contributes to both development of renal insufficiency
and mortality risk in chronic kidney disease. Circulation Research.
2015;116(3):448-455.

Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating
prediction models. Medical Decision Making. 2006;26(6):565-574.

Wang Z, Klipfell E, Bennett BJ, et al. Gut flora metabolism of
phosphatidylcholine promotes cardiovascular disease. Nature.
2011;472(7341):57-63.

## Citing

Citation metadata is in [CITATION.cff](CITATION.cff).

## License

MIT. See [LICENSE](LICENSE).
