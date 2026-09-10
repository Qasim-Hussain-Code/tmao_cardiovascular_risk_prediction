# Statistical analysis plan

This document fixes the analysis before the data are examined. It is written
first so that the choices below cannot be adjusted, consciously or
otherwise, once the results are visible. Any departure from what is written
here should be recorded in the section at the end and reported as a post hoc
decision.

## 1. Objective

To estimate whether fasting plasma trimethylamine N-oxide adds predictive
value to a cardiovascular risk model built from established risk factors,
for the outcome of a major adverse cardiovascular event within three years
of sampling.

The objective is deliberately narrower than asking whether TMAO is
associated with cardiovascular events. Association and predictive usefulness
are different claims, and a marker can be strongly and reproducibly
associated with an outcome while adding almost nothing to a model that
already contains age and renal function. The criteria applied here follow
the American Heart Association statement on evaluating novel risk markers
(Hlatky et al., 2009).

## 2. Population and outcome

Adults with a fasting plasma sample and no prior cardiovascular event at
baseline. The outcome is binary: a major adverse cardiovascular event within
three years, as defined by the source cohort.

The three year horizon is treated as fixed and complete. If a real cohort
has meaningful loss to follow up, a binary outcome model is the wrong tool
and the analysis should move to a time to event framework with
administrative censoring. That is a change of method, not a parameter, and
it would be recorded as a departure.

## 3. Models

Two nested specifications, differing by exactly one term.

**Baseline model.** Age, sex, current smoking, diabetes, body mass index,
systolic blood pressure, total cholesterol, HDL cholesterol, estimated
glomerular filtration rate, and log transformed high sensitivity C reactive
protein.

**Extended model.** The baseline model plus log transformed TMAO.

Renal function is in the baseline model on purpose, and it is the single
most important specification decision here. TMAO is cleared renally, so
impaired kidney function raises circulating TMAO and independently raises
cardiovascular risk (Tang et al., 2015). A model that omitted eGFR would let
TMAO stand in for renal impairment and would overstate its contribution.

The primary estimator is unpenalised logistic regression. Gradient boosted
trees are fitted as a prespecified sensitivity analysis, to check that a
small or absent incremental effect is not an artefact of assuming log TMAO
enters linearly on the log odds scale.

## 4. Handling of the variables

TMAO and C reactive protein are analysed as natural logarithms. Both are
strongly right skewed and both are conventionally modelled that way.

Continuous variables are centred and scaled, so that coefficients are
expressed per standard deviation and are comparable across variables in
different units. The transformation is fitted within each training fold.

Missing values are imputed at the median within each training fold. This is
adequate for a small proportion of missingness that is plausibly unrelated
to the outcome. If a real cohort shows more than about ten per cent
missingness on any analysis variable, or missingness that depends on the
outcome, median imputation understates the uncertainty and the analysis
should move to multiple imputation. That would be recorded as a departure.

## 5. Resampling

Performance is estimated from out of fold predictions under five fold
stratified cross validation repeated five times. Folds are seeded from the
analysis configuration, and every model in a run sees identical folds.

The pairing is a requirement rather than a convenience. The comparison
between the baseline and extended models is paired, and the test used to
compare them assumes both sets of predictions come from the same
participants under the same partitions.

Apparent performance, meaning performance on the data a model was fitted to,
is not reported anywhere. It is optimistic by an amount that grows with the
number of candidate terms, which is precisely the setting under study.

## 6. Outcome measures

**Discrimination.** Area under the receiver operating characteristic curve,
with confidence intervals and a paired comparison between models by the
method of DeLong and colleagues (1988), implemented in the fast form given
by Sun and Xu (2014).

**Calibration.** Calibration intercept and calibration slope from a logistic
recalibration of the out of fold predictions, the Brier score, and the Brier
skill score relative to a model predicting the observed event rate for
everyone. Calibration is reported for every model, not only the primary one,
because a model can discriminate well and still be systematically wrong
about absolute risk, and absolute risk is what a clinician acts on.

**Reclassification.** The integrated discrimination improvement and the
category free net reclassification improvement, each with a percentile
bootstrap interval from one thousand resamples stratified by outcome. The
event and non event components of the reclassification statistic are
reported separately, because a total near zero can conceal a marker that
helps among cases and harms among controls.

These two measures are reported because they are conventional in this
literature, not because they are trusted. Both are sensitive to
miscalibration and both reward a model that is merely more extreme in its
predictions. They are read as descriptive companions to the decision curve.

**Clinical utility.** Net benefit across risk thresholds from one to forty
per cent, against the strategies of treating everyone and treating nobody
(Vickers and Elkin, 2006).

## 7. Primary inferential statement

The primary comparison is the difference in area under the curve between the
extended and baseline logistic models, tested by the paired DeLong test at a
two sided alpha of 0.05.

One comparison is designated primary and the rest are supporting. No
correction for multiplicity is applied, on the grounds that the supporting
measures are described rather than tested. This is a reason to read the
supporting measures as descriptive and not to quote whichever of them looks
most favourable as though it had been the primary question all along.

## 8. What would count as a negative result

A difference in area under the curve whose confidence interval includes
zero, together with a decision curve on which the extended model does not
separate from the baseline model across the plausible threshold range, is a
negative result. It should be reported as such.

Recording this in advance is the point of the section. The most common way a
biomarker study goes wrong is not fabrication but a gradual redefinition of
success once the primary comparison disappoints.

## 9. Departures from this plan

None to date. Any change made after the data were first examined should be
recorded here with its date and reason.
