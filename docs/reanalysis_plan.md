# Reanalysis plan: TMAO and recurrent angina in ST001420

This document is written before any loader, model or estimate exists in the
repository, and it is committed before the code it governs. Its purpose is to
fix the question, the procedures and the interpretation of every possible
outcome while none of them is known. Anything decided after the first estimate
is produced belongs in the departures section at the end, marked as post hoc.

## 1. What changed, and why

The original plan in [analysis_plan.md](analysis_plan.md) asked whether plasma
TMAO adds predictive value to a conventional cardiovascular risk model for
incident major adverse cardiovascular events in a primary prevention setting.
No dataset supporting that question is publicly available. The cohorts that
could answer it, chiefly MESA and the Cardiovascular Health Study, hold their
TMAO measurements behind controlled access.

A search of measured metabolites rather than study titles located one public
dataset that supports a narrower question. The scope change is substantial and
is stated here rather than absorbed quietly into a README.

**The question is no longer** whether TMAO predicts cardiovascular risk.

**The question is now** whether TMAO carries information beyond its own
dietary precursors, in a secondary prevention cohort of post-PCI patients with
stable angina, of Asian ethnicity, uniformly on dual antiplatelet therapy,
predicting recurrent angina within nine months.

That is a different and much narrower claim. It is not a weaker version of the
original question. It is a different question, and nothing learned here
transfers to primary prevention, to hard cardiovascular endpoints, or to
populations unlike this one without further evidence.

The original plan is not deleted. It remains the governing document for the
day a cohort with clinical covariates and adjudicated events becomes available.

## 2. Disclosure of partial unblinding

Before this plan was written, the analyst had already loaded the matrix and
computed univariate summaries of the four pathway metabolites by outcome
group, including a rank sum test on TMAO. That check was performed to
establish that the dataset was usable at all, and its result is known.

This is a real departure from an ideal pre-registration and is recorded here
rather than concealed. Its practical consequence is that the direction of the
univariate TMAO difference is known in advance. The thresholds in section 6
were nevertheless fixed before any model was fitted, any cross validation was
run, or any discrimination, calibration or net benefit statistic was computed.
No multivariable result of any kind has been seen.

## 3. The dataset

Metabolomics Workbench study ST001420, deposited under project PR000974 by
Cui and colleagues at the University of California, San Diego, released 30
July 2020 under CC BY 4.0.

| Property | Value |
| --- | --- |
| Design | Prospective multicentre observational cohort |
| Sampling | Fasting plasma drawn 48 hours after percutaneous coronary intervention |
| Follow-up | Every 30 days to 270 days |
| Participants | 750 hospitalised patients with stable angina after PCI |
| Outcome | Recurrent angina by nine months: 210 events, 540 angina-free |
| Panel | 600 metabolites by targeted LC-MS/MS, SCIEX QTRAP 6500+ |
| Pathway coverage | TMAO, choline, betaine, carnitine, and four related species |

The cohort excludes stenosis above 90 per cent or total occlusion, high
sensitivity C reactive protein above 1.0 mg/dl, serum creatinine above 3.0
mg/dl or dialysis, cancer within five years, prior bypass grafting, iodine
contrast contraindication, and severe hepatic insufficiency.

The renal exclusion is worth noting. Much of the concern about confounding by
renal clearance that motivated putting eGFR in the original baseline model is
partly handled here by design, because patients with substantially impaired
renal function were not enrolled. The corollary is that this dataset cannot be
used to study the TMAO and renal function relationship at all.

The associated publication is Cui et al., Advanced Science, 2021. That paper
reports a multi-metabolite model achieving above 89 per cent accuracy,
sensitivity and specificity across three cohorts. Only the first cohort, the
750 participants analysed here, was deposited. The second discovery cohort of
775 and the validation cohort of 130 were not, so no external validation is
possible from public data, and no one outside the original group has been in a
position to check the reported performance.

## 4. Hard constraint: measurement scale

**The metabolite values in ST001420 are relative peak areas in arbitrary
instrument units. They are not concentrations. They are not micromoles per
litre, nanomoles per litre, or any other physical unit.**

Consequently, and without exception in this analysis:

- No value from this dataset may be compared against any concentration
  threshold, cut point, reference range or tertile boundary drawn from any
  other study.
- No absolute TMAO level reported here may be described as high or low
  relative to any published cohort.
- Only comparisons internal to this dataset are admissible: between
  participants, between outcome groups, and between metabolites measured in
  the same run.

This constraint is enforced in code rather than left to memory. The loader
declares the measurement scale as a typed value, and any operation that would
apply a concentration threshold raises an exception when handed peak area
data. A comment would not survive the months between writing this and someone
reaching for a familiar cut point from the literature. An exception will.

## 5. Hard limitation: acquisition order is confounded with outcome

In the deposited matrix, sample identifiers S1 through S210 are exactly the
210 recurrent angina cases, and S211 through S750 are exactly the 540
angina-free participants. The deposit contains no batch identifier, no run
order, no injection sequence and no acquisition date.

The consequence is structural and cannot be engineered around:

1. If the samples were acquired in the order deposited, any instrument drift
   over the run is perfectly aligned with the outcome, and is indistinguishable
   from biological signal.
2. Because no batch metadata exists, this possibility can neither be tested
   directly nor adjusted for.

This limits what any result from this dataset can mean. In particular, a
finding that a metabolite model discriminates well **cannot** establish that
the discrimination is biological. It establishes only that the discrimination
is reproducible under resampling of these data.

Two consequences follow for interpretation, and they differ by question.

Question 1 compares TMAO against three other metabolites measured on the same
samples in the same run. A drift affecting the whole panel similarly would
largely cancel in that comparison, so question 1 is relatively robust to this
limitation, though not immune, since metabolites differ in their sensitivity
to drift.

Question 2 concerns absolute discrimination of a whole panel model. It is
fully exposed to this limitation. Its result is therefore a statement about
the modelling and validation procedure, not about biology, and it will be
reported as such.

## 6. Questions, procedures, and pre-specified interpretation

All estimates come from out of fold predictions under repeated stratified
cross validation, five folds by five repeats, with folds shared across every
model in a run so that comparisons are paired. Metabolites are analysed as
natural logarithms of peak area. No apparent performance is reported anywhere.

### Question 1, primary. Does TMAO add beyond its own precursors?

Baseline model: log choline, log betaine, log carnitine.
Extended model: baseline plus log TMAO.
Estimator: unpenalised logistic regression.
Primary statistic: difference in out of fold area under the curve, tested by
the paired DeLong test, two sided, alpha 0.05.
Supporting: calibration intercept and slope, Brier skill, integrated
discrimination improvement, category free net reclassification improvement
with stratified bootstrap intervals, and net benefit across thresholds from
0.05 to 0.60, a range chosen to bracket the 28 per cent event rate.

Fixed in advance:

- **Positive.** The confidence interval on the difference in area excludes
  zero, and the extended model's decision curve lies above the baseline
  model's across a substantial part of the threshold range. Conclusion: in
  this cohort, TMAO carries information about recurrent angina beyond its
  dietary precursors.
- **Negative.** The interval includes zero, or the decision curves do not
  separate. Conclusion: TMAO adds nothing beyond its precursors here. This is
  a publishable result and will be stated without softening.
- **Discordant.** Interval excludes zero but the decision curves do not
  separate. Conclusion: a statistically detectable but clinically negligible
  contribution. The decision curve governs the wording.

### Question 2, replication. Does the published discrimination survive honest validation?

Model: L2 penalised logistic regression on all 600 metabolites, with the
penalty selected inside each training fold by nested cross validation so that
no information from a held out fold reaches model selection.
Reported: out of fold accuracy, sensitivity and specificity at the threshold
maximising the Youden index, plus area under the curve, calibration intercept
and slope.

The published claim is above 89 per cent on all three of accuracy,
sensitivity and specificity. Fixed in advance:

- **Reproduced.** Out of fold accuracy, sensitivity and specificity are all at
  or above 0.85, and the calibration slope lies between 0.80 and 1.25.
  Conclusion: the reported discrimination is robust to honest resampling.
  Given section 5, this would still not establish a biological basis.
- **Collapsed.** Any one of the three falls below 0.75, or the calibration
  slope falls below 0.70. Conclusion: the reported performance does not
  survive out of fold evaluation on the deposited cohort. **This is the
  finding, and it will be reported as plainly and as prominently as a positive
  result, and specifically will not be softened because it disagrees with a
  published paper.**
- **Partial.** Values fall between 0.75 and 0.85 with acceptable calibration.
  The numbers are reported with no verdict attached, and the gap from the
  published figures is stated without interpretation.

### Question 3, diagnostic. Is there detectable acquisition order structure?

Motivated by section 5, and specified now rather than reached for if question
2 collapses.

Within each outcome block separately, so that the analysis is blind to the
between group difference, the Spearman correlation between sample index and
log peak area is computed for every metabolite, on complete pairs, requiring at
least 30 non-missing values in the block. Under the null of no order structure
the resulting p values are uniform on the unit interval.

"Materially displaced" is operationalised here rather than judged after the
distribution is seen. Order structure is declared present if **either** the
proportion of metabolites with a within block Spearman p value below 0.05
exceeds 0.10, which is double the null rate, in either block, **or** a
Kolmogorov-Smirnov test of the 600 p values against the uniform distribution
returns p below 0.001 in either block.

- **Order structure present.** Conclusion: acquisition order structure is
  detectable, and any between group difference in this dataset is confounded
  with it to an unknown degree. Every other result in this analysis is then
  reported under that caveat.
- **No order structure detected.** Conclusion: no evidence of drift along the
  deposited ordering.

**Reporting rule for the null branch.** A null result here is the outcome most
likely to be misread. Standing alone beside a clean p value it reads as "the
confound was checked and ruled out", which is false. It was checked under one
unverifiable assumption, namely that the deposited order is the acquisition
order, and the deposit contains nothing that establishes this. Therefore, when
question 3 returns no detected order structure, the limitation from section 5
must be restated in the same breath as the result, in the same sentence and
with equal prominence, never as a footnote or a separate later paragraph. This
is enforced in code: the null verdict string is constructed with the caveat
embedded, and a test asserts that it cannot be emitted without it.

This diagnostic cannot vindicate the dataset. It can only detect a problem or
fail to detect one.

## 7. What will be reported regardless of outcome

Every question above, with its pre-specified verdict applied as written. No
question will be dropped for producing an inconvenient answer, and no
threshold in section 6 will be adjusted after the fact. If a threshold turns
out to have been poorly chosen, the original verdict is reported first and the
reasoning about the threshold follows, labelled post hoc.

The claim scope stated in section 1 will accompany every reported result:
post-PCI secondary prevention, Asian cohort, uniform dual antiplatelet
therapy, symptom-driven endpoint, relative peak areas, single cohort, no
external validation, acquisition order confounded with outcome.

**This scope statement is printed by the reporting code itself**, as the
leading block of the results output and as a header on every results table,
rather than living only in this document where a reader must already know to
look for it. The reason is concrete: a collapsed question 2 verdict sitting
next to an unscoped headline is precisely the thing that gets screenshotted out
of the repository and circulated without its qualifications. The renderer
requires the scope as an argument and verifies that every element listed above
is present in it, raising if any has been dropped. A result cannot be rendered
without its scope.

## 8. Departures from this plan

The partial unblinding described in section 2 is a departure from an ideal
pre-registration and is recorded there.

**Amendments made on 11 September 2026, before any estimate was produced.**
Three changes, all following reviewer comment on the plan as first committed,
and all made while no model had been fitted and no question had been run. They
are amendments rather than departures, and the distinction is that nothing in
the data informed them.

1. Question 3's phrase "materially displaced from the null" was unoperationalised
   in the first version. It is now a stated rule with numbers, fixed before the
   distribution was seen.
2. Question 3 gained an explicit reporting rule binding a null result to the
   section 5 limitation, because a null there is the result most likely to be
   read as exoneration.
3. Section 7 gained the requirement that the scope statement is emitted by the
   reporting code rather than documented here.

**Departures recorded on 11 September 2026, after the first estimates were
produced.** These are post hoc and are labelled as such wherever they appear.

1. **Question 2's partial branch was written for only one of the two ways of
   reaching it.** The plan describes it as values falling between 0.75 and
   0.85. What actually occurred was discrimination far above the reproduction
   floor with a calibration slope of 3.458, far above the acceptable band. The
   verdict label is unchanged, because the decision rule in the code behaved
   exactly as written. The sentence attached to it was factually wrong for this
   case and has been corrected to distinguish the two routes. The thresholds
   themselves were not touched.

2. **Post hoc diagnostics were added**, in `src/tmao_cvd/post_hoc.py`, prompted
   by question 2 returning an out of fold area under the curve of exactly 1.000
   on 750 participants. None was pre-specified. Reporting that figure as a
   performance result and moving on would have been the more serious failure.

3. **The reasoning in section 5 about question 1 was wrong, and this is the
   most important correction here.** Section 5 argued that question 1 was
   "relatively robust" to acquisition artefact because a drift affecting the
   whole panel similarly would largely cancel in a comparison between
   metabolites measured in the same run. The premise fails. The post hoc
   diagnostics show the panel is affected heterogeneously, with univariate
   discrimination ranging from chance to 0.98 across metabolites, a median of
   0.685 where the null is 0.5, and 76.5 per cent of the panel discriminating
   above 0.6. Under heterogeneous drift, one metabolite can outperform another
   simply by being more sensitive to it.

   The consequence is that **question 1's POSITIVE verdict stands as computed
   and pre-specified, but cannot be read as evidence that TMAO carries
   biological information about recurrent angina.** TMAO sits at the 64th
   percentile of the panel with 144 of 405 metabolites ranking above it, and
   betaine, one of its own precursors, discriminates better than it does. An
   incremental contribution measured inside a globally shifted panel is not
   evidence of biology. The verdict is reported unchanged, with this
   qualification attached to it, rather than being quietly revised.

4. **A defect in the panel definition was found and fixed after the first run.**
   `metabolite_columns` excluded the outcome and the sample index but not the
   four derived log transforms, so question 2 ran on 604 columns in which four
   metabolites appeared twice, once raw and once logged. The post hoc
   missingness table reporting 604 rather than 600 metabolites is what exposed
   it. Question 1 was unaffected, since it never uses the panel, and its
   estimates are identical before and after. Question 2's calibration slope
   moved from 3.458 to 2.971 and its verdict was unchanged. A regression test
   now asserts the panel is exactly 600 columns.

Any further change made after this point is to be recorded here with its date
and reason.
