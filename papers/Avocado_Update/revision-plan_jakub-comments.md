# Revision plan — Jakub's comments on *Exploring Equivalence Testing with the Updated TOSTER R Package*

Working notes, drafted 2026-09-04. Line numbers refer to `Avocado_Update.Rmd` at commit `4d24bbc` and will shift once editing starts.

**Summary of my read:** six comments, four of which I think are clearly worth acting on (2, 3, 4, 5), one that deserves a partial/cheap response rather than the restructure it implies (1), and one that is a one-sentence addition with a caveat (6). The heaviest lift by far is #5 (a validation simulation for the ANOVA extension); #4 is the highest value-per-hour.

| # | Topic | Verdict | Effort |
|---|---|---|---|
| 4 | ANOVA section underdeveloped | Accept in full | ~2–3 h writing |
| 3 | Heteroscedasticity "under the null" | Accept, one clause | ~15 min |
| 2 | WMW location-shift vs. no-assumption | Accept, restructure section | ~1–2 h |
| 6 | McNabb & Murayama citation | Accept with caveat | ~30 min |
| 1 | Tempo/genre shifts | Partial — signposting, not restructure | ~1 h |
| 5 | Validate the ANOVA method | Accept — targeted simulation in supplement | ~1 day |

---

## 1. Genre shifts: general stats resource / tool tutorial / new methods

> *"the paper seems to be changing tempo/style between general statistical learning resource, specific tool's tutorial, and presentation of new methods (for ANOVA)"*

**Verdict: real observation, but the fix is not restructuring. Partial accept.**

Pushing back on the implied remedy: this is a software/tutorial paper, and the mixed register is largely intrinsic to the format. The conceptual material — what WMW actually tests, the warning on standardized effect sizes, the estimand footnote in "TOST with t-tests" — is arguably the most valuable content in the paper, and it is what separates this from package documentation. Jakub says as much himself about the effect-size warning. Stripping it out to make a cleaner tutorial would make the paper worse. Splitting the ANOVA method into a separate methods paper would leave the tutorial with a hole in it.

What *is* a genuine problem is that the reader currently has no map. Fixes:

- **Roadmap paragraph at the end of the Introduction** (after the "TOSTER occupies a specific niche" paragraph, ~line 82) that names the three registers explicitly and tells the reader they can skip between them: sections X–Y are conceptual background on choosing an estimand and bounds; sections Y–Z are worked demonstrations organised by design; the ANOVA section describes a generalisation of Campbell & Lakens (2021) that is new to this package.
- **Function-selection table**, early — probably right after "Setting Equivalence Bounds", or as the first thing in "TOST with t-tests". Columns: *research question / estimand*, *design*, *function*, *section*. This is the single highest-return navigational change and doubles as a reference the reader keeps. Rows: mean or mean difference → `t_TOST` / `tsum_TOST` / `simple_htest`; mean with relaxed assumptions → `boot_t_TOST`, `perm_t_TOST`; ratio of geometric means → `log_TOST`; pseudo-median → `wilcox_TOST`; stochastic superiority → `brunner_munzel`; correlation → `z_cor_test`, `boot_cor_test`; SMD → `smd_calc`, `boot_smd_calc`; variance explained / omnibus → `equ_anova`, `equ_ftest`; planning → `power_t_TOST`, `power_z_cor`.
- **Flag the ANOVA novelty in the abstract.** Item (5) currently reads "equivalence testing for ANOVAs", which sounds like a wrapper around existing work. If the generalisation of the non-centrality parameter to factorial and within-subjects designs is new — it is — say so. This also pre-empts comment #4's "most novel but least described" complaint.
- **Optional, larger: promote the estimand framing to a first-class section.** The footnote at ~line 100 is doing a lot of load-bearing conceptual work for a footnote, and it is the thread tying together the rank-based, resampling, log-transform, and ANOVA sections. A short "Choosing an estimand" section sitting beside "Setting Equivalence Bounds" would make the paper read as one argument rather than three.

Explicitly **not** doing: splitting into two papers, or cutting conceptual material.

---

## 2. WMW section jumps between location-shift and no-assumption versions

> *"a more explicit separation via section structure might help people navigate this — in my experience, the awareness of what WMW tests is abysmal"*

**Verdict: accept. He is right, and his read on reader awareness is right.**

The current text (~lines 590–610) interleaves three things in consecutive sentences: what the test does by default, what it does under location shift, and what `wilcox_TOST` puts bounds on. Restructure `### Wilcoxon-Mann-Whitney Tests` into labelled sub-parts:

1. **What WMW tests without extra assumptions.** Stochastic equality; H₀: P(X > Y) + ½P(X = Y) = 0.5. Keep the Fay & Malinovsky (2018) citation here.
2. **What WMW tests under a location-shift assumption.** Only here does the Hodges–Lehmann estimator become a median difference. Keep the counterexamples (`@median_test`) here — their force is precisely that the assumption fails silently.
3. **What `wilcox_TOST` puts bounds on, and what that means in practice.** Bounds go on the pseudo-median. State plainly: *if you are not willing to assert location shift, the pseudo-median is your estimand and your equivalence bound has to be justified on that scale.* Then the existing redirect — `perm_t_test`/`boot_t_test` for means, `brunner_munzel` for stochastic superiority.
4. Worked example, then `#### Checking Assumptions`.

Two supporting changes:

- **One-line estimand statement at the top of each robust-method subsection** (WMW, Brunner–Munzel, bootstrap/permutation, log-TOST), so the parallel structure is visible at a glance. Cheap, and it makes the whole "Robust Methods" section navigable.
- **Consider an ECDF figure for location shift.** Symmetry currently gets a figure (`fig:symmplot`) while location shift — the more consequential assumption in the two-sample case, and the one Jakub is worried readers do not understand — gets a sentence of prose. A two-panel figure (ECDFs as horizontal translations vs. crossing curves) would carry it better. Check `toster_assumption_checks_supplement.qmd` first; there may already be something usable.

---

## 3. "Valid under heteroscedasticity" — specify *under the null*

> *"is it worth specifying the heteroskedasticity is under null? Because outside the null, heteroskedasticity is not an issue?"*

**Verdict: accept. Cheap and correct.**

He is right on the substance. The claim at ~line 645 (`@karch2021`, "should be preferred as the default non-parametric procedure because it is valid under heteroscedasticity") is a Type I error statement: this is the nonparametric Behrens–Fisher problem. WMW's exactness rests on exchangeability, i.e. identical distributions under H₀; when the distributions differ in scale or shape but stochastic equality holds, WMW's error rate is not nominal. Brunner–Munzel tests H₀: p = 0.5 directly and does not require matched distributions.

I would go slightly further than the clause Jakub suggests and name the problem — "nonparametric Behrens–Fisher" is the searchable term, and readers who know the parametric version will get it immediately. Draft replacement:

> `@karch2021` argues it should be preferred as the default non-parametric procedure. Its main advantage over WMW is that the null hypothesis of stochastic equality does not require the two distributions to be identical — the nonparametric Behrens–Fisher problem. WMW's Type I error rate can depart from nominal when the groups are stochastically equal but differ in variance or shape; the Brunner–Munzel test remains valid in that setting and provides a directly interpretable effect size.

Also worth a sentence in the Brunner–Munzel `#### Checking Assumptions` (~line 664). The text correctly says the assumptions are essentially independence plus adequate n for the t-approximation, but the *reason* heteroscedasticity is a non-issue — the estimand is defined from the two marginal distributions, so nothing about it presumes matched shapes — sits a paragraph away from the claim. Tighten.

*(Struck: an earlier version of this note flagged the p̂ direction convention as needing a check, on the mistaken grounds that `mu = c(0.3, 0.7)` was asymmetric. It is symmetric about 0.5, and since swapping X and Y maps p̂ → 1 − p̂, symmetric bounds are invariant to the convention anyway. The convention is also already fixed and documented: the formula method takes X as the first factor level and Y as the second (`R/brunner_munzel.R:1213`), and the print method labels the estimate with the actual level names. No action. The only residual thought — very low priority — is that a reader using genuinely asymmetric bounds, e.g. `mu = c(0.4, 0.8)`, has to know which level is X to place them correctly; if the paper ever shows such an example, name the direction in the surrounding prose. The package's own reporting convention follows `t.test`'s x/y argument style rather than the notation used in the methodological literature, which is a deliberate choice and does not need defending in the paper.)*

---

## 4. ANOVA section: most novel, least described; η² undefined

> *"it's not specified what η^2 is, and overall given that equivalence interpretation for ANOVA is maybe a bit less obvious than some other tests"*

**Verdict: accept in full. Most actionable comment in the email.**

The section (~lines 784–840) is about five paragraphs for the paper's genuinely new methodological contribution, and it assumes the reader already knows partial η², already knows why an omnibus equivalence test is coherent, and already knows how to pick a bound on a variance-explained scale. None of those are safe assumptions.

- **Define partial eta-squared.** η²ₚ = SS_effect / (SS_effect + SS_error) — the proportion of variance attributable to the factor after removing variance attributable to other factors. Say explicitly that this differs from η² (which divides by total SS) and that the two coincide only in a one-way design. Readers pulling values out of `afex` or SPSS output need this.
- **Explain why the test is one-sided, and that it is not a "TOST".** Because η²ₚ ≥ 0 there is no lower bound to test against; the hypotheses are H₀: η²ₚ ≥ Δ vs. H₁: η²ₚ < Δ. This is a non-inferiority / omnibus test, which is exactly why Campbell & Lakens frame it that way. The manuscript never says the two-one-sided-tests logic does not apply here — in a paper built entirely around TOST, that is a real trip hazard.
- **Say what a rejection licenses.** Rejecting means no effect in the design accounts for more than Δ of the (residual-adjusted) variance. That is a statement about the factor as a whole, not about any particular pairwise contrast, and it does not license "the groups are equivalent" in the mean-difference sense. Conversely, failing to reject the equivalence null is not evidence of a large effect. Two or three sentences.
- **Bound selection on the η²ₚ scale, which is harder than on the raw scale.** Two problems worth naming: (a) partial η² is not comparable across designs, because the denominator depends on which other factors are in the model and, for within-subjects designs, on which error term the effect is tested against — so a bound borrowed from a published study in a different design is not portable; (b) the Δ = 0.35 used in *both* worked examples is very large and is there only to make the demo produce a non-degenerate result. The paper has a global "these are illustrative" disclaimer at the end of "Setting Equivalence Bounds", but 0.35 appearing in two consecutive examples reads as a default and needs a local repeat. Consider giving the conversion to Cohen's *f* (f² = η²ₚ/(1−η²ₚ)) as a bridge for readers with *f*-scale intuition, and pointing at the bound-selection section of the `the_ftestTOSTER` vignette.
- **State the generalisation in the manuscript, not only on the website.** Put both formulas inline: Campbell & Lakens' one-way λ = NΔ/(1−Δ), and the generalised λ_eq = [Δ/(1−Δ)]·(df₁ + df₂ + 1). Then note that for a one-way design df₁ + df₂ + 1 = (J−1) + (N−J) + 1 = N, so the generalisation reduces *exactly* to the published result. That reduction is free, it is a real (if partial) validation, and it answers both #4 and part of #5. It is currently nowhere in the manuscript.
- **Walk the `equ_anova` output columns.** The existing sentence lists `p.null`, `pes`, `eqb`, `p.equ` in running prose; a short definition list would be easier to use. The within-subjects example is a good place to note that each effect gets its own error term and therefore its own denominator for η²ₚ.

---

## 5. Was the original ANOVA method validated? Add validation here.

> *"I'm not sure how much the original work on equivalence testing for ANOVA validated the method, but it feels like something that might work nicely if added here."*

**Verdict: accept, with a correction to the premise. This is the big one.**

Correcting the premise: Campbell & Lakens (2021) *did* validate their method — the paper includes a simulation study of Type I error and power for the one-way ANOVA and multivariable regression R² cases, plus a comparison against a Bayesian alternative. So the published method is not unvalidated, and the manuscript should say so (one clause, currently absent, worth adding regardless).

What is **not** validated is my generalisation. The λ_eq = [Δ/(1−Δ)]·(df₁ + df₂ + 1) substitution is what extends the method to factorial and within-subjects designs, and it reduces to the published formula only in the one-way case. Honest framing: the base method is validated, the extension is not. A reviewer is likely to find this too, so better to close it now.

Proposed scope — a simulation supplement (probably its own file rather than folded into `toster_assumption_checks_supplement.qmd`, since it is a different kind of content):

- **Designs:** (a) one-way between-subjects, as a check that the implementation reproduces Campbell & Lakens; (b) 2×2 between-subjects factorial — main effects and interaction; (c) 2×2 fully within-subjects, which is the case the generalisation is least obviously right for, since the ncp for a repeated-measures F depends on the correlation among measures and it is not self-evident that df₁ + df₂ + 1 carries the right information there; (d) mixed design, if (c) behaves.
- **Conditions:** true η²ₚ below, at, and above the bound — the boundary is the least favourable case under H₀ and gives the Type I error; a few n per cell (small / medium / large); for the within-subjects cases, a range of correlations among repeated measures; Δ spanning a plausible range (say 0.02 / 0.06 / 0.14) rather than the 0.35 used for demos.
- **Outcomes:** empirical rejection rate at the boundary (should be ≤ α, and ideally close to α rather than badly conservative), power curves, and a direct comparison of the generalised p-value against the Campbell & Lakens formula in the one-way case — should be numerically *identical*, not merely close.
- **Deliverable:** one figure (rejection rate vs. true η²ₚ, panelled by design) plus a short table, referenced from the ANOVA section with a single summary sentence in the manuscript. Detail stays in the supplement.

Add as regression tests in `tests/testthat/` once the simulation exists: exact equality of `equ_ftest` against the closed-form one-way result, and a fixed-seed boundary Type I error check for a factorial case.

**Risk to flag:** if the within-subjects case turns out to be miscalibrated, that is a package bug rather than a paper problem, and it changes the schedule. Run condition (c) first as a smoke test before committing to the full grid.

---

## 6. Cite the "multilevel modelling is often unnecessary" paper

> *"it might be worth mentioning this paper showing that quite often multilevel-handling-statistics is unnecessary and simple averaging can do a decent job"*

**Verdict: accept, but with a caveat and in a specific place.**

The reference is McNabb & Murayama (2021), *Unnecessary reliance on multilevel modelling to analyse nested data in neuroscience: When a traditional summary-statistics approach suffices*, **Current Research in Neurobiology** 2, 100024, doi:10.1016/j.crneur.2021.100024. Note for the .bib: the ScienceDirect PII Jakub sent resolves to *Current Research in Neurobiology*, not *NeuroImage: Reports*.

It fits naturally in the **Conclusions** limitations paragraph (~line 894), which currently notes that TOSTER does not handle multilevel/mixed-effects designs and points to `emmeans`/`marginaleffects`. Adding this citation turns a bare limitation into useful advice: for nested designs of the common "several trials per participant, roughly balanced" kind, averaging within participant and running a simple test on the per-participant summaries is often defensible — which brings a fair number of nominally multilevel designs back inside TOSTER's scope.

**Caveat to include rather than citing it flat:** there is a published commentary — Bloom, Thieu & Bolger (2022), *Current Research in Neurobiology* — arguing that the summary-statistics approach breaks down under unequal cluster sizes or unequal cluster variances (it over-weights small clusters), and that it cannot recover between-cluster heterogeneity even when it gets the fixed effect right. That is directly relevant to *equivalence* testing specifically: an equivalence conclusion drawn from a summary-statistics analysis of unbalanced clusters could be driven by mis-weighting rather than by a genuinely small effect. So cite McNabb & Murayama, state the balanced-design condition, and cite the commentary in the same sentence. One or two sentences total — do not open the debate in a package paper.

Minor scope pushback: the McNabb & Murayama argument is framed around neuroscience trial-level data, and this paper's audience is broader. Keep the phrasing conditional ("where clusters are of similar size and the question concerns the average effect") rather than implying multilevel modelling is generally avoidable.

---

## Consolidated to-do

**Manuscript (`Avocado_Update.Rmd`)**
- [ ] Intro: roadmap paragraph; flag the ANOVA generalisation as novel (also in the abstract).
- [ ] Add function-selection table.
- [ ] (Optional) Promote the estimand footnote to a short section.
- [ ] Restructure WMW subsection: no-assumption / location-shift / what-gets-bounded.
- [ ] One-line estimand statement at the head of each robust-method subsection.
- [ ] (Optional) ECDF location-shift figure alongside the symmetry plot.
- [ ] Brunner–Munzel: qualify "valid under heteroscedasticity" as a null/Type I error statement; name the nonparametric Behrens–Fisher problem.
- [ ] ANOVA: define partial η²; distinguish from η²; explain the one-sided non-inferiority framing and why it is not a TOST; state what a rejection does and does not license.
- [ ] ANOVA: bound-selection discussion; flag Δ = 0.35 as illustrative; add the f² conversion.
- [ ] ANOVA: both ncp formulas inline, plus the one-way reduction df₁ + df₂ + 1 = N.
- [ ] ANOVA: note that Campbell & Lakens ran a validating simulation; reference the new simulation for the extension.
- [ ] ANOVA: definition list for `equ_anova` output columns.
- [ ] Conclusions: McNabb & Murayama with the balance condition, plus the Bloom et al. commentary.

**Supplement / package**
- [ ] New simulation supplement for the factorial and within-subjects extension (run the within-subjects smoke test first).
- [ ] `tests/testthat/`: exact-equality test against the one-way closed form; boundary Type I error check for a factorial design.

**Bibliography (`interactcadsample.bib`)**
- [ ] `mcnabb2021`, `bloom2022` entries.

**Reply to Jakub**
- [ ] Note that Campbell & Lakens did validate the one-way case, and that his comment correctly identifies the extension as the unvalidated part — this is the most useful thing in his email.
- [ ] Say what I am not doing on #1 (no restructure) and why.
