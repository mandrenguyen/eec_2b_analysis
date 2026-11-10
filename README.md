# eec_2b_analysis

EEC (energy-energy correlator) between the two b-jets in double-b events, PPRef2024.

**`workflow/` is the live code.** Everything at the top level and in `Run3/`,
`TemplateFit_Run3/`, `Unfolding/` is an older copy — read those, do not run them.

Run everything on LLR (`ssh llruicms01.in2p3.fr`), from
`/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow`. On the Mac the same tree is the
sshfs mount: `/Users/zoezaidan/llruicms01/...` == `/home/llr/cms/zaidan/...`. See `CLAUDE.md`
for the mount and the LCG/Condor environment.

This is the only readme: a quick start first, then the full reference behind every step. The
archived `fnb` / `fb` / `lnfb` sections are in `workflow/readme_archive_fraction_observables.md`.

## Setup

```bash
source setup_roounfold_env.sh          # LCG 106a + the private RooUnfold build
root -l -b -q -e ".L create_files_for_template_fit.cpp++"   # compile the .so once
```

## ⚠️ Recent changes — read before running anything

### 2026-10-01 — jet energy corrections, whole chain re-run

- **JEC in step 1.** `create_files_for_template_fit.cpp` can recompute `jtpt` from `rawpt` with the 2024 pp-reference JEC (L1FastJet + L2Relative on MC and data, + L2Residual on data), exactly as Matt's `*_JEC_Prompt24HIpp` productions do — checked jet by jet, identical. Switched by the new last argument `apply_jec` (default off) and by `APPLY_JEC=true` in both run scripts, which then require an `OUT_TAG` containing `jec`. The correction files come from the `Run3HeavyIonJetCalibrations` submodule (`git submodule update --init` if it is empty).
- **Everything re-produced with it** as `_upartv2rhojec` / `_upartv2rhojec_B` (drivers: `condor/run_rhojec_production.sh`, `condor/run_rhojec_analysis.sh` on `/data_CMS`). `jecInputs()` in `result_paths.h` is **on**, so every macro reads the JEC results; set it false to read the uncorrected `_upartv2rho` ones, which are untouched.

### 2026-09-28 — both observables re-produced on the `_rho` inputs

- **Input era.** Marc's chunks are now `merged_block_*_UParTV2_rho.root` (the old ones are gone). EEC(ΔR) and dN/dB are both re-produced on them, so they come from the same MC. The switch is `rhoInputs()` in `result_paths.h`, and everything carries `eraTag()` = `_upartv2rho`: the productions (`OUT_TAG=upartv2rho` for ΔR, `upartv2rho_B` for B), the template-fit directories and the result folders. The pre-`_rho` results are untouched; set `rhoInputs()` to false to read them. ΔR data is now our own `_upartv2rho` production (no longer Afnan's).
- **Run scripts.** Both read `INPUT_SUFFIX` (default `_rho`). `SKIP_COMPILE=true ./run_agg_ntuple_chunks.sh` reuses the built `.so`, so several productions (and the data jobs) can run at once. Data job scripts go in one directory per production; submit with the `condor_submit condor.submit -append 'Executable = …'` line the script prints.
- **Same band on both.** In the rho era the ΔR band books exactly B's four sources: MC template modelling, detector response, tracking efficiency and the light-jet mistag pair. UParT SF and EEC-weight modelling are out of the ΔR band. The UParT SF *correction* is still applied to the ΔR nominal.
- **Labels.** Every plot says **CMS Work In Progress**, MC-only plots included. The pre-fit first-look plots (1–4) now carry the header too.

### 2026-09-25 — fewer files

The two readmes (this file and `workflow/readme_workflow.md`) are now this one file.

**Four macros were folded into existing files**, to keep the number of files down. The code is unchanged; only the command line is different (load the host file, then
call the function):

| was | now lives in | run it as |
|---|---|---|
| `plot_z_first_look.C` | `template_fit.cpp` | `root -l -b -q -e '.L template_fit.cpp+' -e 'plot_z_first_look("qcd","pythia","upartv2_B")'` |
| `plot_systematics_summary.C` | `apply_weights_and_systematics.C` | `root -l -b -q -e '.L apply_weights_and_systematics.C' -e 'plot_systematics_summary("both")'` (and `plot_final_generator_band`) |
| `momentum_balance_mc_study.C` | `plot_purity_efficiency_response.cpp` | `root -l -b -q -e '.L plot_purity_efficiency_response.cpp' -e 'momentum_balance_mc_study("both","pythia")'` |
| `compare_root_contents.C` | `Help_Functions.h` | `root -l -b -q -e '#include "Help_Functions.h"' -e 'compare_root_contents("a.root","b.root")'` |

The systematics talk (`workflow/talk_systematics/`) is kept locally only, outside git.

### 2026-09-24 — the measurement is EEC(ΔR) and dN/dB, each with the full band

The B band now books every systematic the ΔR one has, except the two that cannot exist for B
(UParT SF, EEC-weight modelling). See *What is in each band* and *The full B band*. ⚠️ Its
inputs are **not produced yet**, and the upstream input chunks changed under us: read the
blocked note in *The full B band* before launching anything.

`apply_unfolding_2d.C`'s `TF_GENERATOR` default is now `"pythia"` (it was `"herwig"`), so the
nominal call really is the nominal.

The archived `fnb` / `fb` / `lnfb` sections moved to `readme_archive_fraction_observables.md`.

### 2026-09-22 — two changes

**The momentum balance is called `B`, not `z`.** The rename went all the way down: the
`OBSERVABLE` argument is `"B"`, the histograms carry the suffix `_B`, the results paths carry
`_B`, and the production tag is `_upartv2_B`. Nothing translates between the two spellings,
so **the whole chain has to be re-run for the balance** — step 1 into `OUT_TAG=upartv2_B`,
then the fit, the unfolding and the band. Pointing this code at a pre-rename production
(`_zfirst`, `_3obs`, `_lnfblin`, `_fb`) fails on a missing `h3D_data_B`, which is the right
way round but is not a fallback.

The `_z` files and result folders already on disk are **not** deleted, renamed or overwritten.
They stay reproducible from the code that made them. Every path in this readme has been
rewritten to the **new** `_B` spelling, so a path quoted below is where a re-run will put its
output, not necessarily where the old run's output is sitting today.

`plot_z_first_look()` keeps its name, but its default observable is now `"B"`.

**Only ΔR (EEC) and `B` are measured.** `fnb`, `fb` and `lnfb` are commented out; their
sections are kept for the record in `readme_archive_fraction_observables.md`.

## Quick start

| Step | Script | Produces |
|---|---|---|
| 1. MC chunks | `run_agg_ntuple_chunks.sh` | per-block `RMatrix_*`, `*MCGEN*`, `AggBHadronNtuple_*` |
| 2. Data chunks | `make_hardprobes_condor_scripts.sh` + `condor_submit` | per-block data templates |
| 3. Template fit | `template_fit.cpp` | signal fraction per bin, nominal and varied fit files |
| 4. Unfolding | `apply_unfolding_2d.C` | unfolded EEC + closure/scan plots |
| 5. Systematics | `apply_weights_and_systematics.C` | final result + band, `data_vs_gen_*`, `systematics_curves_*` |

Steps 1 and 2 both run the same macro, `create_files_for_template_fit.cpp`. Step 3 fits what
they wrote, step 4 unfolds it, and step 5 reads one step-4 output per systematic variation.

### 1. MC ntuples — `run_agg_ntuple_chunks.sh`

- Runs `create_files_for_template_fit.cpp` over the MC chunks, one background ROOT job per block, producing the response matrices, gen histograms and aggregated B-hadron ntuples.
- Edit `SAMPLE` (`qcd`|`bjet`) and `GENERATOR` (`pythia`|`herwig`) at the top, or override them on the command line.
- Takes ~1 h and prints nothing while it runs; never run two copies at once.

```bash
./run_agg_ntuple_chunks.sh
SAMPLE=qcd GENERATOR=herwig ./run_agg_ntuple_chunks.sh
TRACK_EFF_UNC=true ./run_agg_ntuple_chunks.sh     # tracking-efficiency systematic
```

### 2. Data ntuples — `make_hardprobes_condor_scripts.sh`

- Writes the 50 per-job Condor scripts that run the same macro over HardProbes data; usually nothing to edit.
- `condor.submit` is the submit file they are launched with.

```bash
./make_hardprobes_condor_scripts.sh
condor_submit condor.submit
condor_q <cluster> -af:h JobStatus RemoteHost    # 1 idle, 2 running, 5 held
```

### 3. Template fit — `template_fit.cpp`

- Fits the 2B signal and the effective background to data and writes the signal fraction per ΔR bin, plus the nominal and varied (`var0B_2`, `var0B_0`) fit files.

```bash
root -l -b -q 'template_fit.cpp("both","pythia")'
```

### 4. Unfolding — `apply_unfolding_2d.C`

- Sums the per-block files of steps 1–2, applies the signal fraction and the UParT SF to data, unfolds and applies the purity and efficiency corrections.
- **Matrix inversion is the default** (`unfoldBayes = false`), for both observables: it has no regularisation, so there is no per-observable iteration count to justify (Bayesian needed 7 for ΔR and 21–26 for B). Pass `true` for the Bayesian result.
- Nothing to edit — every knob is an argument: `(SAMPLE, UNFOLDING_GENERATOR, test_mode, unfoldBayes, scan_niter, TF_GENERATOR, TRACK_EFF_UNC, TF_VARIATION, SFUPART_VARIATION, EEC_WEIGHT_OFF, OBSERVABLE, NITER)`; `test_mode` is `0` full-MC closure, `1` split test, `2` data; `OBSERVABLE` is `dr` or `B` (the momentum balance). The three b/non-b fraction observables (`fnb`, `fb`, `lnfb`) were commented out on 2026-09-22; see *Recent changes* above.
- Any observable other than `dr` forces `SFUPART_VARIATION = "off"` (the UParT SF is binned in reco ΔR and has no equivalent) and `EEC_WEIGHT_OFF = true` (it is a yield: dN/dB). Both rules live in `result_paths.h`, one function each.
- One results folder per argument combination, under `/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/`. **Matrix-inversion results go in their own subdirectory**, `results/matrix_inversion/unfolding_*`, so a whole MI production — every variation, both observables, plots and logs — is one directory to sync or compare. Bayesian keeps the flat layout it always had, so nothing that exists today moved.

```bash
root -l -b -q 'apply_unfolding_2d.C("both","pythia")'                  # nominal
root -l -b -q 'apply_unfolding_2d.C("both","herwig")'                  # detector response
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",true)'               # tracking eff
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_2")'    # mistag up
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_0")'    # mistag down
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","jpcalib_hf")'  # UParT SF
root -l -b -q 'apply_unfolding_2d.C("both","pythia",0)'                # closure test
# the momentum balance (dN/dB): nominal, then its two mistag variations
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","off",true,"B")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_2","off",true,"B")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_0","off",true,"B")'
```

### 5. Systematics and final plots — `apply_weights_and_systematics.C`

- Reads one step-4 result per booked variation, builds the systematic band, and writes the final ROOT file and the two plots you show; a missing variation is reported and skipped.
- Outputs go into the nominal folder of that observable: `final_*_with_systematics.root`, `data_vs_gen_*`, `systematics_curves_*`, `sfupart_before_after_*`.
- Run it **once per observable**: `OBSERVABLE` is the last argument and picks both the results it reads and the variation list. ΔR gets the full set (template fit, detector response, tracking, mistag pair, UParT SF stat + calibration); B gets the mistag pair, because the Herwig and track-drop MC have no B axis in them yet — the other entries are in the file, commented, with what each needs.
- To add a systematic: one entry in the `variations` vector, and its tag in `result_paths.h`.

```bash
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'   # EEC(ΔR)
root -l -b -q 'apply_weights_and_systematics.C("both","pythia",2,false,"h_data_fully_corrected_2D","pythia",2,false,"B")'   # dN/dB
```

## Other macros

| File | What it does |
|---|---|
| `plot_purity_efficiency_response.cpp` | Plots the response matrix, purity and efficiency from the step-1 files. Also holds `momentum_balance_mc_study()`, the purity/stability study that chose the B binning. |
| `template_fit.cpp` → `plot_z_first_look()` | The pre-fit step for any observable (last argument): flavour templates, m_2B separation per bin, purity/efficiency, data vs MC, and the raw yields. |
| `apply_weights_and_systematics.C` → `plot_systematics_summary()`, `plot_final_generator_band()` | Plots every Herwig variation as a relative shift against the Pythia nominal (needs both nominals from step 5). |
| `plotNice_UParT_roc.C` / `plotNice_bdt_roc.C` | ROC curves for the UParT b-tagger and for the BDT. |
| `merge_RMatrixTH2D_qcd.C` / `_bjet.C` / `_bjet_qcd_merged.C` | Merge per-block response matrices into one file (step 4 sums in memory, so this is only for inspection). |
| `Help_Functions.h` → `compare_root_contents()` | Bin-by-bin comparison of two ROOT files, to prove a rerun reproduces exactly. |
| `sync_plots_to_eos.sh` | Copies plots from `/data_CMS` to the CERNBox website (needs `kinit zzaidanc@CERN.CH`). |

## Headers

| File | What it holds |
|---|---|
| `result_paths.h` | The one definition of every result path and variation tag; add a variation here and both macros follow. |
| `binning_histos_small.h` | The ΔR, jet-pT and mass binning. |
| `central_selections.h` | The jet and event selections. |
| `tracking_efficiency_syst.h` | The hashed 3% track drop used by the tracking systematic. |
| `tTree.h` | The input ntuple branches. |
| `Help_Functions.h`, `Draw_EEC.h` | Plotting helpers and the EEC drawing used by the template fit. |
| `rootlogon.C` | Sends ACLiC build products to `/data_CMS`; loaded automatically by `root` in this folder. |

---

# Reference, step by step

Everything below is the detail behind the quick start.

## 1. MC — `run_agg_ntuple_chunks.sh`

**Change (usually only these two):**

| Line | Knob | Values |
|---|---|---|
| 13 | `SAMPLE` | `qcd` \| `bjet` |
| 14 | `GENERATOR` | `pythia` \| `herwig` |

Everything else follows: input dir, output dir, chunk filename pattern, the macro's sample
tag, and the list of blocks (taken from what is on disk — 10 qcd/pythia, 9 bjet/pythia,
8 qcd/herwig, 9 bjet/herwig).

Occasionally:

| Line | Knob | Note |
|---|---|---|
| 18 | `INPUT_TAG` | `UParTV2` \| `negTagFix` \| `negTag`. Also picks `BTAG_WP`. Only Pythia8 QCD has the `negTag*` chunks. |
| 38 | `OUT_TAG` | suffix on every output file; change it to avoid overwriting an earlier run |
| 56 | `TRACK_EFF_UNC` | `false` nominal, `true` the tracking-efficiency systematic — see below |

**Run:**

```bash
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && ./run_agg_ntuple_chunks.sh'
```

It prints the resolved config first — check it before walking away:

```
sample bjet herwig (UParTV2), b-tag WP 0.712, sample tag bjet
input  /data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/bJetHerwig/UParTV2_chunks
output /data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/bJetHerwig/agg_ntuple_chunks
found 9 chunks
```

It compiles the macro once, then runs one `nice`d ROOT job per block in parallel and waits.
Outputs land in `<OUT_BASE>/block_NNNN/`, logs in `<OUT_BASE>/logs/block_NNNN_<OUT_TAG>.log`.

**It takes ~50–60 min and looks hung the whole time.** The macro prints
`Processing events [0, N)` once and the staging files only grow at autoflush boundaries. To
check it is alive, watch the counters rise:

```bash
for p in $(pgrep -u zaidan -f "root.exe.*create_files_for_template_fit"); do
  awk '{print $14+$15}' /proc/$p/stat; grep ^read_bytes /proc/$p/io
done
```

Do **not** kill and restart (no resume — each block starts from event 0), and **never run
two copies at once**: it does `rm -rf` on the staging dirs and wipes the shared ACLiC build.

### Tracking-efficiency systematic (`TRK_DROP`)

Track reconstruction in the jet core is mismodelled in simulation. To propagate that,
set `TRACK_EFF_UNC=true`: 3% of the **reconstructed** tracks in MC are thrown away during the
B reconstruction (`makeSvtxs_withBDT`, before any vertex is built), so the SVs, m_2B, ΔR and the EEC weight are all built from a jet
that is genuinely missing tracks. The gen-level b hadrons (`refTrk*`) are untouched — the
variation has to move reco while truth stays put, which is what makes the response matrix
and the corrections change. Data is never varied (the macro guards on `isMC`).

**It is reproducible by construction.** The number that decides a track's fate is not
drawn from a stream, it is *hashed* from `(TRK_SEED, entry number, track index)`. So the
same track always gets the same verdict, no matter how many blocks the sample is split
into, in what order they run, or how many times you rerun. A seeded `std::mt19937`
consumed as a stream would **not** give this: it only reproduces if the tracks are
consumed in exactly the same order by exactly one process, which the parallel block jobs
violate. See `tracking_efficiency_syst.h` for the argument in full.

Verified on `block_0000` of Pythia8 QCD: two runs of the same varied job give bit-identical
histogram contents (46 + 13 histograms, max |Δ| = 0), and `[0,50k) + [50k,100k)` summed
equals a single `[0,100k)` run up to 1e-16 rounding.

Output files are tagged `_trkdrop030`, so a variation run cannot overwrite the nominal and
both can sit in the same `block_NNNN/` directory. Nominal filenames are unchanged.

**The full procedure:**

```bash
# 1. produce the varied MC — TRACK_EFF_UNC=true in run_agg_ntuple_chunks.sh, once per
#    SAMPLE/GENERATOR you use in the nominal (qcd and bjet if you unfold with "both")
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && ./run_agg_ntuple_chunks.sh'

# 2. unfold the SAME data with it (last argument), which lands in its own results folder
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",true)'

# 3. the band: the "tracking_eff" entry in apply_weights_and_systematics.C picks the
#    varied result up automatically and adds |nominal - varied| to the quadrature sum
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'
```

The signal fraction **is** re-derived on the varied templates: with `TRACK_EFF_UNC=true`,
`apply_unfolding_2d.C` reads the fit from `TemplateFits_<sample>_pythia_trkdrop030…/`, so
`template_fit.cpp(SAMPLE, GENERATOR, true, …)` has to have run on the hadd'ed track-drop
MCGEN files first. That refit is what dominates this systematic on EEC(ΔR). Holding the
signal fraction at its nominal value instead is one line: drop `track_eff_unc` from the
`templateFitFile()` call in `apply_unfolding_2d.C`.

### Light-jet mistagging systematic (`TF_VARIATION`)

The **0B template** is the jets that passed the b tag with no gen b hadron in them
(`jtNbHad == 0`) — mistagged light and charm jets. It is not a free component of the
template fit: the fit has two PDFs, 2B (signal) and an *effective* background, and 0B is
folded into that background at the ratio MC predicts,

```
c' = 0B / (0B + 1B)          effective bkg = (1-c') x 1B_shape + c' x 0B_shape
```

so nothing in the data constrains the mistag rate. It reaches the measurement only through
the **shape** of the effective background, and what moves is the fitted signal fraction.

The systematic rebuilds that shape with the 0B admixture scaled by ×2 and ×0 — a
deliberately conservative ±100% on the mistag rate, which brackets the measured
light-flavour mistag scale factors (typically 20–50% at a tight WP) by a wide margin and
needs no external SF input. `template_fit.cpp` already runs all of them: every fit
directory holds `nominal_`, `var0B_2_` and `var0B_0_` copies of
`Run3_TemplateFits_histos_3d_80_inf.root`, each with its own `h_sig_fraction_fit`.

**So there is no MC to produce and no refit — only two more unfoldings:**

```bash
# 1. the two varied signal fractions, each into its own results folder
#    (_tf0Bx2 / _tf0Bx0)
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_2")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"var0B_0")'

# 2. the band: the "mistag_0B_up" / "mistag_0B_down" entries pick them up automatically
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'
```

The two entries share `Variation::group = "mistag_0B"`, so the total takes their **per-bin
envelope once** instead of adding both in quadrature — they are two sides of one
uncertainty, not two sources. Any other up/down pair can be combined the same way by giving
it a shared group.

`TF_VARIATION` only makes sense with `test_mode = 2`: the signal fraction is applied to data
alone, so the macro refuses the closure modes rather than write a folder that reads as a
result and is not one.

**Size, measured on `both`/`pythia` (Pythia templates, UParT WP 0.712):** the shift on the
signal fraction is **≤ 2.4 %, and under 1 % in most bins** — smaller than the fit's own
statistical error everywhere except the first two ΔR bins. It is the smallest of the
systematics booked so far. The largest shifts sit at the two ends: the lowest-ΔR bins, where
the sample is background-dominated (signal fraction 0.12–0.14) so the background shape
matters most, and the highest-ΔR bins, where `c'` itself is largest (up to 0.18).

There is a guard in `template_fit.cpp`: `eff_bkg1B = 1 - eff_bkg0B`, so `c' > 0.5` would
make the ×2 variation build a background PDF with **negative** 1B weight. The normalisation
check further down would still pass — the two weights sum to 1 by construction — and
`RooHistPdf` would silently clip the negative bins. It cannot happen at the current WP
(`2c' ≤ 0.37`), but a looser b tag or coarser ΔR bins can reach it, so the value is now
clipped with a warning that says the variation is truncated in that bin.

⚠️ `c'` is computed from the **dijet (qcd) integrals only** — `h_bb`/`h_b`/`h_nob`, never
the `_bjet` ones — while the shapes being mixed are qcd+bjet when `also_bjet`. That is
deliberate (qcd is the sample whose flavour composition matches the data; the bjet sample is
there for template statistics and would bias `c'` downwards if it entered the ratio), but it
is worth knowing when reading the numbers.

## 2. Data — `make_hardprobes_condor_scripts.sh`

**Usually change nothing.** It only ever runs HardProbes data, so `SAMPLE_TAG=data`
(line 19) is fixed. Occasionally line 15 `INPUT_TAG` or line 22 `OUT_TAG`.

⚠️ The `BTAG_WP` case here is duplicated from `run_agg_ntuple_chunks.sh` — if you change a
working point, change it in **both** files or data and MC come out at different WPs.

**Run — two steps, and the `.so` must exist first:**

```bash
# 0. compile once (run_agg_ntuple_chunks.sh also does this); jobs exit 3 if it is missing
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && \
  source setup_roounfold_env.sh && root -l -b -q -e ".L create_files_for_template_fit.cpp++"'

# 1. generate the 50 job scripts
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && ./make_hardprobes_condor_scripts.sh'

# 2. submit
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && condor_submit condor.submit'

# 3. monitor  (JobStatus: 1 idle, 2 running, 5 held; gone from the queue = finished)
ssh llruicms01.in2p3.fr 'condor_q <cluster> -af:h JobStatus RemoteHost'
```

Read the job logs locally through the mount rather than over SSH.

## 4. Unfolding — `apply_unfolding_2d.C`

**Nothing to edit — everything is a macro argument:**

| # | Argument | Values | Default |
|---|---|---|---|
| 1 | `SAMPLE` | `qcd` \| `bjet` \| `both` | `qcd` |
| 2 | `GENERATOR` | `pythia` \| `herwig` | `pythia` |
| 3 | `test_mode` | `0` full-MC closure, `1` split test, `2` data | `2` |
| 4 | `unfoldBayes` | `true` Bayesian, `false` matrix inversion | **`false`** |
| 5 | `scan_niter` | `true` scan niter 1..100, pick the best by refolding GoF | `false` |
| 11 | `OBSERVABLE` | `dr` \| `B` (momentum balance) | `dr` |

### The two MC tests — `test_mode` 0 and 1

Both use **MC only**: no data file and no template fit, so an observable can be validated
before its data production exists. Run them for every observable before quoting a result:

```bash
for m in 0 1; do for o in dr B; do
  root -l -b -q "apply_unfolding_2d.C(\"both\",\"pythia\",$m,false,false,\"pythia\",false,\"nominal\",\"nominal\",false,\"$o\")"
done; done
```

`SFUPART_VARIATION` and `EEC_WEIGHT_OFF` are passed as `"nominal"`/`false` above on purpose —
the macro forces each observable's own convention from `nominalSfupart()` /
`nominalEecWeightOff()` and says so in the log, so the caller never has to know them.

⚠️ **The refolding GoF is not the test here.** Matrix inversion reproduces its input exactly,
so refolding returns `chi2/ndf = 0.000, p = 1.000` in *both* modes for *every* observable —
including a split test that failed. What each mode actually tells you:

| mode | what it checks | how to read it |
|---|---|---|
| `0` full closure | response, corrections and truth from the same events | `h_data_fully_corrected_2D / h_mc_true_2D` must be **1.0000 in every bin**. Anything else means the response and the corrections were not built from the same events — a bug in the chain, not physics. |
| `1` split test | corrections + response from the odd half, pseudodata + truth from the even half | the real test, and it is allowed to fail. Compare the same two histograms bin by bin. |

⚠️ **The split test compares against the gen of the *pseudodata*, not the gen that filled the
matrix.** The two halves are cut on `ient % 2` in `create_files_for_template_fit.cpp`:

| half | what it fills | histogram |
|---|---|---|
| even (`ient % 2 == 0`) | the pseudodata, **and its own gen truth** | `h3D_pseudodata_bb`, `h_pseudodata_truth_tf` |
| odd (`ient % 2 == 1`) | the response matrix and the purity/efficiency ratios | `response_tf_pseudo_full`, `h_full_pseudo_*_tf` |

`apply_unfolding_2d.C` unfolds the even half with the odd half's response and compares the
result to `h_pseudodata_truth_tf` — the even half's truth. Comparing instead against
`h_full_pseudo_efficiency_denominator_tf` (the odd half's gen) would be testing the unfolding
against the prior it was built from: it would close better and mean nothing. **Both halves
are ~half the sample, so a half-swap would not show up as a scale error** — it has to be read
off the histogram name, which is why the legend now spells out which half every curve is.

The bottomline plot is drawn (`unfolding_plot_*_bottomline_test_eec_2D.pdf`) but the numbers
are not logged — read them off `h_data_fully_corrected_2D` and `h_mc_true_2D` in the output
`histos_*_after_unfolding_2D.root`.

**Results, `both`/`pythia`, matrix inversion, 2026-09-21.** Full closure is exact (ratio
1.0000 in every bin) for all three. Split test:

| observable | bins | χ²/ndf vs truth | p | worst bin |
|---|---|---|---|---|
| `dr` | 9 | 11.89/9 = 1.32 | 0.220 | bin 8, ratio 0.954 at 2.3% stat |
| `B` | 5 | 1.19/5 = 0.24 | 0.946 | bin 3, ratio 1.040 at 4.2% stat |
| `lnfb` | 5 | 5.57/5 = 1.11 | 0.351 | bins 1–2, ratio 0.82/0.87 at 31%/12% stat — see the archive |
| `fb` | 5 | 12.09/5 = 2.42 | **0.034** | bin 1, ratio **0.524** at 28% stat — **FAILS**, see the archive |

That χ² counts only the unfolded statistical error. The two halves are independent, so adding
the truth error would lower it — read these as upper bounds.

### Bayesian iteration counts

`NITER` is the macro's **12th** argument and is **not written to any log**, so a Bayesian run's
regularisation strength cannot be recovered from its output. These were measured by
`scan_niter` on data (refolding GoF, full 2D reco space) on 2026-09-21:

| observable | niter (full 2D) | niter (pT bin 2 only) |
|---|---|---|
| `dr` | 7 | — |
| `B` | 25 | 20 |
| `lnfb` | 37 | 29 |
| `fb` | 21 | 19 |

Two numbers per observable, each tuned on its own scan, is exactly why matrix inversion is the
nominal unfolding: a band whose variations are unfolded at a regularisation strength of their
own is not a band. Use these only to reproduce the Bayesian cross-check.

Bayesian closure and split tests were run for all three observables at these counts on
2026-09-21 — they had never existed before that date for any observable.

### A curated copy of the whole chain

`~/Documents/plots/final_results/` on the Mac holds the four steps for all three observables,
both unfolding methods, with a README carrying the numbers above. Rebuild it by re-copying
from the results area; it is a view, not a source.

A wrong `SAMPLE`/`GENERATOR` string is rejected with a message instead of silently falling
back to the default.

### Why matrix inversion

Matrix inversion is the unfolding of **both** observables. It has no regularisation, so
there is no iteration count to justify per observable — the Bayesian runs needed `niter = 7`
for ΔR and 21–26 for B, two numbers each tuned on its own refolding scan, and a band whose
variations are unfolded at a regularisation strength of their own is not a band. `scan_niter`
is refused with `unfoldBayes = false` rather than silently ignored, and the log says which
unfolding ran instead of printing an `niter` that means nothing.

Pass `unfoldBayes = true` to get the Bayesian result back.

### Where each production lives

```
results/                                  <- Bayesian, the flat layout it always had
  unfolding_both_pythia_upartv2/
  unfolding_both_pythia_trkdrop030_upartv2/
  ...
results/matrix_inversion/                 <- one MI production, complete
  unfolding_both_pythia_upartv2/
  unfolding_both_pythia_B_sfupartoff_noeecw_upartv2/
  ...
```

`methodDir()` in `result_paths.h` is the whole rule: `""` for Bayesian, `matrix_inversion/`
for inversion. It is a **directory** rather than another name tag because it is not a
variation of the measurement — the same variation unfolded two ways is the same systematic,
and interleaving two complete productions in one parent is what makes a results area
unreadable. Every plot, log and `final_*` file follows the result it belongs to, so a whole
production can be synced to the website, diffed or deleted as one thing.

Bayesian paths are untouched, including the ones already used in talks. Files still carry
`_bayesian` / `_MI` in their names as well, so a file is self-describing even once it has
been copied out of its folder.

⚠️ `resultFolder()` takes `unfoldBayes` as its **third** argument, before the defaulted ones,
so a call written against the old signature fails to compile instead of quietly building a
Bayesian path for an inversion run.

**What it costs, measured** (EEC(ΔR), `both`/`pythia`, data, 2026-09-17). Inversion is
unregularised, so it neither pulls towards the prior nor damps the statistical fluctuation,
and both effects land in the first two bins, where the migration is largest:

| bin | Bayesian (niter 7) | matrix inversion |
|---|---|---|
| 1 | 1.755 ± 0.061 stat, syst 22% | 2.168 ± 0.115 stat, syst 27% |
| 2 | 2.776 ± 0.040 stat, syst 8% | 2.233 ± 0.104 stat, syst 23% |
| 3–9 | central values, syst 1.6–9.6% | central values within 6%, syst 3.7–16.3% |

Bins 3–9 keep their central values (−3% to +6%) and carry a somewhat larger band. Bins
1–2 move in the central value too, so read those two as inversion telling you how much of
the Bayesian result there was coming from the regularisation.

### The fraction observables `fnb`, `fb`, `lnfb` — archived

Disabled 2026-09-22. Their sections (definitions, binning studies, results and the `fb`
split-test failure) are in `readme_archive_fraction_observables.md`.

### The `B` observable

`OBSERVABLE = "B"` measures the momentum balance between the two aggregated B hadrons. Two
things about it are properties of the observable, not of a run, so `apply_unfolding_2d.C`
forces them and says so in the log — both rules are single functions in `result_paths.h`:

- `nominalSfupart("B")` → `"off"`: the UParT SF is measured in reco ΔR bins and has no B
  equivalent, so every B path carries `_sfupartoff`. A **known gap**, not a silent skip: the
  ΔR result carries a correction (0.933 → 1.149 across its range) that the B result does not.
- `nominalEecWeightOff("B")` → `true`: B is a yield, so the EEC weight is off and the result
  is dN/dB. Forced because a weighted-B production does exist, and a B result unfolded with
  the weight on would land in a folder nothing reads.

**Run** (on LLR, after `source setup_roounfold_env.sh`):

```bash
root -l -b -q 'apply_unfolding_2d.C("qcd","pythia")'      # data, matrix inversion
root -l -b -q 'apply_unfolding_2d.C("both","herwig")'     # qcd+bjet Herwig
root -l -b -q 'apply_unfolding_2d.C("bjet","herwig",0)'   # full-MC closure
```

In one line from the Mac (note the nested quoting):

```bash
ssh llruicms01.in2p3.fr 'cd /home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow && \
  source setup_roounfold_env.sh && root -l -b -q "apply_unfolding_2d.C(\"both\",\"herwig\")"'
```

Results go to `/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/matrix_inversion/unfolding_<SAMPLE>_<GENERATOR><tags>_upartv2/`
(no `matrix_inversion/` for a Bayesian run) — one directory per flag combination, named by
`resultFolder()` in `result_paths.h` — with a timestamped `.log` per run (stdout is
redirected there, so the terminal stays quiet).

The macro reads the **per-block** files of step 1 directly and sums them in memory; there is
no merge step to run first. `SAMPLE=both` appends the bjet blocks to the qcd ones.

### How "both" is summed

- **count histograms** and the **`RooUnfoldResponse`** are added (`TH1::Add`,
  `RooUnfoldResponse::Add`) — the same operation as `merge_RMatrixTH2D_bjet_qcd_merged.C`;
- **purity/efficiency ratios are never added.** They are recomputed from the summed counts.
  Adding the N stored ratios gives N× the true value (17 Herwig blocks peak at ~16.9).

Sanity check: with `test_mode=0` the closure must be exact — refolding `chi2/ndf = 0.000,
p = 1.000`. If it is not, the response and the corrections are not from the same events.

Physics caveat: `both` adds a b-enriched sample onto inclusive QCD, which double-counts
b-jets unless the sample weights already handle the overlap. Comparing `qcd` vs `both`
shows how much it moves.

The **template fit** file (signal fraction, needed only for `test_mode=2`) is not produced
here — it is read from `results/TemplateFit_Run3/TemplateFits_<both|qcd>_<TF_GENERATOR><tags>_upartv2/`,
written by `template_fit.cpp` (`templateFitFile()` builds the path). `TF_GENERATOR` picks it,
not `GENERATOR`: the two are separate systematics. `SAMPLE=both` takes the `both` fit,
anything else the `qcd` fit.

---

## 5. Systematics and the final plots — `apply_weights_and_systematics.C`

Reads one step-4 result per variation, builds the band, and writes the two plots you
actually show. **Nothing to edit for a normal run** — the variations are the
`std::vector<Variation> variations` list near the top of the macro, and each entry names the
unfolding run it needs.

Run it **once per observable**. `OBSERVABLE` is the last argument and picks both the results
it reads and the variation list; `unfoldBayes` must match the runs (`false` = matrix
inversion, the default on both sides).

```bash
root -l -b -q 'apply_weights_and_systematics.C("both","pythia")'                             # EEC(dr)
root -l -b -q 'apply_weights_and_systematics.C("both","pythia",2,false,"h_data_fully_corrected_2D","pythia",2,true)'        # dN/dr, the yield run
root -l -b -q 'apply_weights_and_systematics.C("both","pythia",2,false,"h_data_fully_corrected_2D","pythia",2,false,"B")'   # dN/dB
```

### What is in each band

**The two measurements are EEC(ΔR) and dN/dB** (decided 2026-09-24). The dN/dΔR yield run
still works but is not a deliverable; its column is kept for reference.

| Source | EEC(ΔR) | dN/dB | dN/dΔR (not a deliverable) |
|---|---|---|---|
| MC template modelling (Herwig template fit) | ✔ | ✔ | — |
| Detector response (Herwig unfolding MC) | ✔ | ✔ | — |
| Tracking efficiency (−3% tracks, MC and refit) | ✔ | ✔ | ✔ |
| Light-jet mistag (0B ×2 / ×0) | ✔ | ✔ | — |
| UParT SF stat ±1σ | ✔ | n/a — no UParT SF exists in B | — |
| UParT SF JP HF swap, qq rate ±25% | ✔ | n/a | — |
| EEC-weight modelling | ✔ | n/a — weight off | n/a — weight off |

The B band is the ΔR band minus the two rows that cannot exist for B, booked with the same
entry names, groups and `sysLabel`s so the two can be compared source by source. The UParT SF
rows are a **known gap** rather than a missing run: the SF is binned in reco ΔR and has no B
equivalent, so the ΔR result carries a correction the B result does not. How to produce every
B input is in *The full B band* below.

A variation whose result file is missing is reported with the exact `apply_unfolding_2d.C`
call that would produce it, and skipped — the nominal alone still works, you just get a
smaller band.

**Outputs**, all into the *nominal* result folder
(`unfolding_both_pythia_upartv2/`):

| File | What it is |
|---|---|
| `data_vs_gen_<label>.{pdf,png}` | the money plot: unfolded data with stat bars and the systematic box, against particle-level MC |
| `sfupart_before_after_<label>.{pdf,png}` | the UParT SF on its own: the result before and after it, and their ratio. Only produced while the display-only `sfupart_off` entry is booked |
| `systematics_curves_<label>.{pdf,png}` | **per-source breakdown** — one curve per variation, and below it each source's fractional pull against the total band |
| `final_<label>_with_systematics.root` | `h_nominal`, `h_syst_total_up`/`_down`, `h_syst_total`, `h_result_stat_syst`, `g_syst_band` (asymmetric), `h_syst_*`/`h_delta_*`/`h_var_*` per source, and `h_gen_*` |

### Where the results are, and which of them are results

The results area gets one folder per flag combination, and **most of them are inputs to the
band, not results**. `apply_weights_and_systematics.C` prints a `FINAL OUTPUTS` block at the
end of every run naming the deliverables, so this never has to be worked out by hand.

All under `/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/results/matrix_inversion/` for the
current (matrix-inversion) production, and directly under `results/` for the Bayesian one —
the folder names below are the same in both:

| Folder | What it is |
|---|---|
| `unfolding_both_pythia_upartv2/` | **the result.** Nominal run, and where all four deliverables are written |
| `unfolding_both_herwig_upartv2/` | input — "Detector response" |
| `unfolding_both_pythia_tfherwig_upartv2/` | input — "MC template modeling" |
| `unfolding_both_pythia_trkdrop030_upartv2/` | input — "Tracking efficiency" |
| `unfolding_both_pythia_tf0Bx2/x0_upartv2/` | input — "Light jet mistagging" pair |
| `unfolding_both_pythia_sfupartjphf/qqup/qqdn_upartv2/` | input — UParT SF systematics (currently commented out of the band) |
| `unfolding_both_pythia_sfupartoff_upartv2/` | comparison — the before-SF run, display only |

The four deliverables, all in the nominal folder:

| File | What it is |
|---|---|
| `final_<label>_with_systematics.root` | the result: central values, stat, σ⁺/σ⁻, per-source histograms |
| `data_vs_gen_<label>.{pdf,png}` | the measurement vs particle-level MC |
| `systematics_curves_<label>.{pdf,png}` | per-source breakdown |
| `sfupart_before_after_<label>.{pdf,png}` | the UParT SF on its own |

⚠️ Three folders predate the current naming and are **not** produced by this code any more:
`unfolding_upartv2/`, `unfolding_both_pythia_respherwig_upartv2/`,
`unfolding_both_herwig_tfherwig_upartv2/`. Nothing reads them. Left in place rather than
deleted, but do not mistake them for current results.

⚠️ `/data_CMS` is **not** on the sshfs mount, so none of this opens from the Mac. Use
`sync_plots_to_eos.sh` or `scp`.

### Paths and variation tags live in `result_paths.h`

Every result path and every variation tag has **one** definition, in `result_paths.h`, which
both `apply_unfolding_2d.C` (which writes them) and `apply_weights_and_systematics.C` (which
reads them) include.

This is not tidiness — it is a bug fix. The two macros used to carry private copies of the
same mapping, and on 2026-09-11 a new variation (`"off"`) was added to one copy only. The
reader's copy did not recognise it, returned `""`, and `""` **is** the nominal tag — so it
loaded the nominal result and compared it against itself. Nothing errored; the column just
read `-0`, which looks like "the correction does nothing" rather than "wrong file".

Three guards now, in increasing order of how early they catch it:

1. **Validation up front.** `apply_weights_and_systematics.C` checks every booked
   variation's sample/generator/tfVariation/sfupartVariation against `isKnown*()` before
   opening anything, and refuses to run naming the bad entry. It also rejects two variations
   sharing a name, which would silently overwrite each other's histograms.
2. **Unknown ⇒ poison, never nominal.** An unrecognised value returns `kUnknownTag()`
   (`_UNKNOWN_VARIATION`), not `""`, so the path cannot alias onto the nominal — it simply
   fails to open, with an error naming the bad value.
3. **The writer checks the reader's path.** After each unfolding, `apply_unfolding_2d.C`
   rebuilds the path `apply_weights_and_systematics.C` will look for and warns if it does not
   exist — so writer/reader drift surfaces at write time, not days later.

**To add a variation, add it in one place: `result_paths.h`.** Add the value to the relevant
`isKnown*()` and its tag to the matching `*VarTag()`, and both macros follow.

### The UParT b-tag efficiency scale factor

The b-tag efficiency SF for the UParT selection on jets with **at least two reconstructed
secondary vertices**, from `lifetime_jp_sfb_agg_dr_wp712_eec_systematics.root` in this
folder, binned in reco R_BB over the same 9 ΔR bins as the analysis.

**It is applied inside `apply_unfolding_2d.C`, to the reco-level data, before unfolding** —
right after the signal fraction multiplies `h_data_after_fit`. That is the point of it: the
SF is measured in *reco* R_BB, so it belongs on a distribution that is still in reco R_BB,
and the response matrix then carries it through the migration. Applying it to the unfolded
result instead would put a reco-binned correction on particle-level bins and silently ignore
that migration.

It is **divided** in: the chain corrects data with the efficiency measured in MC, and the
true efficiency is `ε_data = SF · ε_MC`, so `data/ε_MC` has to be divided by SF again. Data
only — it never touches MC, and the closure modes are refused. It runs 0.933 at low ΔR to
1.149 at bin 8, so it is a real reshaping, not a wash.

The central SF is a **correction, not an option**: every data run applies it, and an untagged
results path means "with the nominal SF", not "without one".

⚠️ **The KEYS and the histogram names differ in that file.** Key `h_SFb_dr_central` holds a
histogram whose internal name is `h_SFb_dr_wp712_eec_central`. `TFile::Get()` takes the
**key**, so the internal name returns null. `upartSFHist()` uses key names.

`SFUPART_VARIATION` (9th argument) picks which SF to divide in, and tags the results folder.
The name is namespaced to UParT on purpose — more scale factors are coming, and each should
get its own `sfupart*`-style plumbing rather than sharing one generic `SF_` knob:

| `SFUPART_VARIATION` | Key | Tag | Treatment in the band |
|---|---|---|---|
| `nominal` | `h_SFb_dr_central` | — | the correction itself |
| `jpcalib_hf` | `h_SFb_dr_jpcalib_hf` | `_sfupartjphf` | one alternative calibration → one-sided, symmetrised |
| `qqrate_up` | `h_SFb_dr_qqrate_up` | `_sfupartqqup` | grouped pair with ↓, envelope once, stays asymmetric |
| `qqrate_down` | `h_SFb_dr_qqrate_down` | `_sfupartqqdn` | ditto |
| `statup` | `h_SFb_dr_central` **+ its bin error** | `_sfupartstatup` | grouped pair with ↓ — the calibration's statistical precision |
| `statdn` | `h_SFb_dr_central` **− its bin error** | `_sfupartstatdn` | ditto |
| `off` | — | `_sfupartoff` | no SF at all; display-only comparison, never in the band |

The `statup`/`statdn` shift is **coherent** across ΔR — every bin moves together. That is
deliberate: the analysis's last two ΔR bins share one calibration bin, so their errors really
are 100% correlated, and a coherent shift gets that right while staying conservative
elsewhere.

```bash
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","jpcalib_hf")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","qqrate_up")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","qqrate_down")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","statup")'
root -l -b -q 'apply_unfolding_2d.C("both","pythia",2,false,false,"pythia",false,"nominal","statdn")'
```

`SFUPART_VARIATION = "off"` applies **no** SF at all. It is not a systematic — it exists so
the same macro produces the before-SF result (tag `_sfupartoff`) for a like-for-like
comparison. It is booked in `apply_weights_and_systematics.C` as a **display-only** entry
(`Variation::inBand = false`): drawn and printed, but contributing nothing to σ⁺/σ⁻ and
nothing to the tag, because it is a correction, not an uncertainty. It drives two things —
the "Before UParT SF" curve on `systematics_curves_*`, and the dedicated
`sfupart_before_after_*` plot (the two curves and their ratio, nothing else). Comment that
one entry out and both disappear.

`sfupartVarTag()` lives once, in `result_paths.h`, and an unrecognised value returns
`_UNKNOWN_VARIATION` rather than `""` — see *Paths and variation tags* above for why that
matters.

The SF's **statistical** error is **not** folded into the data's error. Applying the SF is a
pure rescale — content and the data's statistical error multiplied by the same 1/SF, so the
relative error is unchanged and the error bar stays a *data* error bar — and the SF's own
precision is booked as the `statup`/`statdn` systematic instead. Doing both would count it
twice.

`h_SFb_dr_syst_up`/`_syst_down` are the quadrature of the two calibration systematics and are
**not** used: booking `jpcalib_hf` and `qqrate_up`/`down` separately gives the same total
while keeping each visible as its own column and curve.

⚠️ Because the SF now lives inside the unfolding, **every** unfolding run changes when it
changes — the nominal and all six systematic variations have to be re-run together, not just
the SF ones.

### MC-derived corrections: bin-by-bin rescaling

Purity (before unfolding) and the three post-unfolding corrections — reconstruction
efficiency, 2SV+b-tag efficiency, EEC weight — are applied by `applyCorrection()` as a plain
**bin-by-bin rescale**: content and error multiplied by the *same* factor, so the data's
**relative** error is unchanged.

Deliberately **not** `TH2::Multiply`/`Divide`. Those propagate the correction's own error
into the result, which folds the MC's limited statistics into the *data's* statistical error.
The convention here is to treat every MC-derived correction as exact. The helper also runs
over underflow/overflow, checks the binning explicitly instead of failing silently, and
reports any bin where a zero correction emptied non-empty data.

⚠️ **The MC statistical uncertainty is therefore not counted anywhere.** It is out of the
stat error and no systematic is booked for it — a deliberate choice (2026-09-11), not an
oversight. It is *not* negligible: removing it cut the statistical error by 13–22%, which
backs out to an MC contribution roughly as large as the data's own. Revisit if the error
budget is ever questioned.

⚠️ Corrections applied at **fill time** — the EEC weight and sample weights inside
`create_files_for_template_fit.cpp` — cannot get this treatment. They are baked into the
histogram contents with their fluctuations already folded in, so there is no separable factor
left to rescale, and their MC statistics stay mixed into the data's error.

### The band is asymmetric

σ⁺ and σ⁻ are computed and stored separately, because a genuine up/down variation need not
move the result by the same amount each way. Which side a source feeds depends on whether it
is *really* two-sided:

- **Grouped entries** (a shared `Variation::group` — an up/down pair) are treated **signed**.
  The member that pushes the result up sets σ⁺, the one that pushes it down sets σ⁻, and the
  source contributes its largest shift on each side once.
- **Ungrouped entries** have one alternative (a Pythia↔Herwig swap, the 3% track drop). One
  alternative tells you the *size* of the difference, not its direction — nothing says the
  truth lies on the far side of the nominal — so these are **symmetrised**: |Δ| into both
  sides.

With only the mistag pair genuinely two-sided today, the two bands come out equal to ~3
significant figures. **That is expected** — the machinery is there for the asymmetric sources
still to come, and a new one is asymmetric the moment it is booked as a group.

If both members of a pair happen to push the same way in a bin, that source would contribute
nothing to one side, which claims zero uncertainty there. Too aggressive, so such a group
falls back to its envelope on both sides and prints a `NOTE:` naming itself.

`h_syst_total` is still written as the **larger** of the two sides, for anything that wants
one number per bin (`h_result_stat_syst`'s bin errors, and any script predating the split).

⚠️ **Two sign conventions, deliberately.** The printed table shows `variation − nominal`, the
shift of the *result*, so its sign matches which of σ⁺/σ⁻ the source feeds. The stored
`h_delta_*` keep the original `nominal − variation`, because `plot_systematics_summary()`
reads them and flips the sign itself. Same numbers, opposite sign.

### The per-source breakdown plot

`systematics_curves_*` is where you read **what each systematic does on its own**. The
bottom pad is the content: `variation / nominal − 1` for every source, one line each, over
the shaded total band — so a source's line reaching the edge of the band means it dominates
that bin, and a line sitting on zero means it does nothing there. The top pad is the same
variations as absolute EEC curves, which mostly shows that they all sit on top of the
nominal; it is a reference, not the message.

The pads are split 55/45 rather than the usual 70/30 for exactly that reason — the pulls are
a few percent and need the room. Up/down pairs appear as **two separate lines** here (not
the envelope that enters the total), which is what you want for debugging: it shows whether
a pair is genuinely two-sided or whether one side does all the work.

Per-bin numbers for the same thing are in the printed table — one signed column per source —
and in `final_*_with_systematics.root` as `h_delta_<name>` (signed shift), `h_syst_<name>`
(the symmetrised `|shift|` that enters the band) and `h_var_<name>` (that variation's own
EEC curve).

### The data-vs-gen plot

Top pad: red = the measurement (unfolded data, stat error bars, systematic as a shaded box
per point), blue = particle-level MC of the nominal generator, purple dashed = the other
generator when its unfolding exists. Bottom pad: **MC / Data**, with the data's own
uncertainty drawn around 1 — the systematic as the wide shaded box, the statistical error as
the narrow hatched one, so the two are separable by eye.

Upper left is a **"Systematic uncertainties" tag** listing the sources in the band. It is
built from the `sysLabel` field of the entries that *actually entered* the band, deduplicated
so an up/down pair is listed once — a source whose result file was missing and got skipped
cannot be claimed on the plot. Leave `sysLabel` empty to keep an entry out of the tag while
it still counts in the total. The second particle-level curve is off by default; flip
`show_alt_gen` to draw the other generator next to the nominal one.

The gen curve is `h_mc_true_2D`, which in data mode is `hgenjet_2b_all`: the particle-level
reference the corrected data is meant to land on. It is only re-normalised to unit area,
never weighted — the scale factors in `Weights` correct *data* for detector effects and have
no meaning at particle level.

⚠️ **`h_mc_true_2D`'s `Write()` was disabled until 2026-09-10**, so any result unfolded
before then does not carry it and the plot is skipped with a message. Re-run
`apply_unfolding_2d.C` for that configuration to get it.

### Not booked yet

The band is deliberately these five for now. Known gaps, recorded so they are not
rediscovered from scratch:

| Missing | Status |
|---|---|
| Template **fit range** | `fitRange0to8` is already produced by `template_fit.cpp` on every run (it is inside the `ivar < 4` loop) and consumed nowhere. Now that `TF_VARIATION` exists it is two more unfoldings to book — the cheapest of these by far. `FITRANGE_0_7` is in the enum but the loop stops at 4, so it is not even produced. |
| Unfolding **regularisation** | Closed, by removing the regularisation rather than varying it: matrix inversion is the default for both observables, so there is no iteration count in the result at all. The price is visible in the first two ΔR bins — see below. |
| **JES / JER** | Absent from the whole chain. Moves jets across the pT bin edges, so it does not fully cancel in a unit-area shape. |
| b-tag **efficiency** SF | Only the light-jet *mistag* is covered. `Weights::sf_file` is empty, so no SF is applied at all. A dR-flat efficiency cancels in the normalisation; a dR-dependent one does not. |

⚠️ **Double counting between two that ARE booked.** "Detector response" (`unfolding_model`)
swaps the *entire* unfolding MC to the other generator, which includes the EEC weight
correction — so it already contains "MC modeling of EEC weight". Booking both counts that
difference twice. Read it as the deliberately conservative option, or drop
`unfolding_model` and book the corrections individually.

### Adding a systematic

One entry in the `variations` vector; everything else — the quadrature total, the printed
table, both plots, the written histograms — iterates that vector. The fields are the
`apply_unfolding_2d.C` arguments of the run that produced it:

```cpp
{ "name", SAMPLE, generator, tfGenerator, "legend text", trackEffUnc, tfVariation, group, "Tag name" }
```

`group` is the only non-obvious one. Leave it `""` for an independent source and its
symmetrised `|nominal - variation|` enters the quadrature sum on its own. Give **two**
entries the *same* non-empty group when they are the two sides of one uncertainty (an
up/down pair): the total then takes their per-bin **envelope once** instead of adding both,
which would inflate that source by √2. `mistag_0B_up` / `mistag_0B_down` are the worked
example.

---

## Macro arguments (`create_files_for_template_fit.cpp`)

```cpp
create_files_for_template_fit(RunN, pT_low, etaCut, n, btag, isMC, btagWP,
                              makeTemplates, createRmatrix, makeAggNtuple,
                              ev_first, ev_last, inputFileOverride,
                              outputFolderOverride, sampleTag,
                              track_eff_unc, trkSeed)
```

`track_eff_unc` defaults to `false` (nominal, a no-op) and `trkSeed` to 20260908; the run
script sets them from `TRACK_EFF_UNC` / `TRK_SEED`. The 3% itself is
`TrkEffSyst::kDropFraction`.

There is no `dataType` any more: **`isMC`** is data vs MC, and **`sampleTag`**
(`data`/`bjet`/`qcd`, set by the run scripts) is what names the output histograms
`..._3D_<sampleTag>_f`. `pT_high` is also gone — the upper pT edge lives in
`jtpt_binsVector` in `binning_histos_small.h`.

---

## A second observable: the momentum balance `B`

`B = pT_lead / (pT1 + pT2)` ∈ [0.5, 1), built from the **same two aggregated B hadrons**
whose pT product is the EEC weight. `B = 0.5` is a balanced pair, `B → 1` a very asymmetric
one. It is a **parallel** measured axis, ΔR-style — templates `TH3D(m_2B, B, pT)`, fit
sliced in `(B, pT)`, response 2D `(B, pT)` — not an extra dimension on the ΔR measurement.

**Status (2026-09-24): every step of the chain handles it, including the full band** — see
*The full B band* at the end of this section.

### How the second axis works

The second histogram axis used to *be* ΔR — the name was written into every histogram name,
loop bound and binning lookup by hand. It is now a parameter, `ObsDef` in `observables.h`,
and `create_files_for_template_fit.cpp` fills **one histogram set per observable in a single
pass over the data**. Adding `B` therefore costs no extra run time: the same job that makes
the ΔR templates and response makes the `B` ones.

| | ΔR | B |
|---|---|---|
| `ObsDef::suffix` | `""` | `"_B"` |
| histogram names | `h3D_bb`, `response_tf_full`, … | `h3D_bb_B`, `response_tf_full_B`, … |
| bins | 9, `dr_binsVector` | 4, `B_binsVector` |
| reco-level cut | `dr > 0.005` | none |

**ΔR's suffix is empty on purpose**, so every ΔR object keeps the exact name it has always
had and nothing downstream had to change. That is verified, not assumed — see below.

`observables.h` also holds the single definition of `B` (`MomBalance::value`), shared with
`momentum_balance_mc_study()`. It is deliberately **not** a function pointer inside `ObsDef`:
ΔR comes from `tTree::calc_dr`, which truncates through `Float_t` twice, and reimplementing
that to satisfy a uniform interface would move ΔR in its last bits for no reason. The caller
computes the value; `ObsDef` owns the binning, the overflow fold and the naming.

⚠️ `obsList()` and `obsValues()` in `create_files_for_template_fit.cpp` are a **matched
pair** — the Nth value must belong to the Nth observable. Adding a third observable means
editing both; a length mismatch aborts the job at startup rather than filling one
observable's values into another's histograms.

⚠️ The `dr > 0.005` cut gates the **whole template block**, so a jet it drops is dropped from
`B` too. Deliberate: the two measurements must share a reco domain, or the signal fraction
fitted on one does not apply to the other.

### The ΔR result is unchanged — verified

The refactor was A/B tested on 100k events of `block_0000`, QCD/Pythia, before and after:

```
RMatrix …_qcd_f.root   46 histograms, 0 differ, max |delta| = 0  ->  IDENTICAL
…_qcd_fMCGEN.root      13 histograms, 0 differ, max |delta| = 0  ->  IDENTICAL
all 6 RooUnfoldResponse objects (Mresponse element by element)   ->  IDENTICAL
```

`compare_root_contents()` only compares `TH1`s, so the six response objects were compared
separately — that is the check that would have caught a response filled on a different axis.

Cross-check that `B` is filled on *exactly* the same jets: `h3D_bb` and `h3D_bb_B` integrate
to the same number (to 3e-16, one ulp — the same fills summed over different bins).

To repeat it after a change here:

```bash
# baseline: stash the current source, run 100k events, keep the output
root -l -b -q -e "gSystem->Load(\"$WORKFLOW/create_files_for_template_fit_cpp.so\"); \
  create_files_for_template_fit(3,80,2,1,true,true,0.712,true,true,false,0,100000, \
  \"<a merged_block_*.root>\",\"/tmp/ab_before\",\"qcd\",false,20260908,false)"
# then the same into /tmp/ab_after with the new source, and:
root -l -b -q -e '#include "Help_Functions.h"' -e 'compare_root_contents("/tmp/ab_before/<f>","/tmp/ab_after/<f>")'
```

### The binning study

`momentum_balance_mc_study()` (in `plot_purity_efficiency_response.cpp`) decided the binning **before** the reprocessing was spent:

```bash
root -l -b -q -e '.L plot_purity_efficiency_response.cpp' -e 'momentum_balance_mc_study("both","pythia")'
```

**It reprocesses nothing.** `AggBHadronNtuple` already stores `recoPt1/recoPt2` and
`genPt1/genPt2`, so reco and gen `B` are already on disk for MC. It reads the per-block
nominal ntuples directly (same discovery rule as `aggChunkFiles()`), applies the
response-matrix gates verbatim (`nRecoAgg==2 && passRecoKin && passBtag && recoDr>0.005`,
`nGenAgg>=2 && passGenKin`) and the same `w_reco = weight * recoEec` weight, so its
purity/stability mean what they will mean in the real chain. Outputs land in
`results/momentum_balance_study_<sample>_<generator>_upartv2/`.

⚠️ The purity/stability it prints are **diagonal fractions** of the migration matrix
(`N(reco in bin AND gen in bin) / N(reco in bin)`), *not* the `purity`/`efficiency`
histograms the chain unfolds with — those count how many jets pass reco **and** gen
selection, not how many stay in their bin. Different question, same word.

### What it found (both/pythia, 3.1M jets in the response)

The macro prints ΔR's own numbers in the analysis's 9 bins as the **benchmark**, because a
purity of 0.4 is disqualifying if ΔR scores 0.85 and unremarkable if ΔR scores 0.37:

| | worst bin purity | worst bin stability | typical bin |
|---|---|---|---|
| ΔR, the 9 analysis bins | 0.354 | 0.373 | ~0.60 |
| `B`, 4 equal-occupancy bins | **0.417** | **0.404** | ~0.50 |
| `B`, 5 equal-occupancy bins | 0.360 | 0.342 | ~0.40 |

So on the equal-occupancy scan, **4 bins clears the ΔR bar and 5 does not**.

**The binning actually in use is 5 uniform bins of 0.1** (chosen 2026-09-15), not the
equal-occupancy set — a round-number binning was preferred:

```
B_binsVector = { 0.500, 0.600, 0.700, 0.800, 0.900, 1.000 }
```

Scored on `both`, that is four good bins and one weak one:

```
0.5-0.6   purity 0.610   stability 0.573    33% of entries
0.6-0.7          0.462             0.455    29%
0.7-0.8          0.425             0.452    21%
0.8-0.9          0.387             0.446    12%
0.9-1.0          0.251             0.216     5%   <-- below the ΔR bar
ΔR worst         0.354             0.373
```

⚠️ **The 0.9–1.0 bin reconstructs more migrated-in jets than its own** (purity 0.251 against
ΔR's worst of 0.354). It holds only 5% of the entries, so it costs little elsewhere, but do
not read that bin as a measurement without checking it against the `test_mode=0` closure.
Cause: `B_reco` is pulled toward the bulk (~0.67) from both ends, and at `B_gen ~ 0.99` the
bias reaches −0.175, larger than the 0.100 resolution itself.

⚠️ **Changing the B binning forces a full re-production.** The binning is baked into the
histograms at fill time, and these edges do not align with the previous equal-occupancy set,
so nothing can be rebinned after the fact. Budget ~1h per MC sample plus the condor data run.

**The resolution is flat but the bias is not.** `RMS(B_reco − B_gen) = 0.100` across the
whole range, and it barely moves with `B` (0.074 → 0.105). What moves is the **mean**:

```
B_gen ≈ 0.51  ->  mean(B_reco - B_gen) = +0.069
B_gen ≈ 0.70  ->                          0.000
B_gen ≈ 0.99  ->                         -0.175
```

Reco `B` is pulled towards the bulk (~0.67) from both ends — ordinary regression to the
mean from a 0.1-wide smearing kernel on a bounded variable, and at high `B` it is **larger
than the resolution itself**. That is why the top bin has to be wide: 0.76–1.00 is about
two bias-lengths, and anything narrower is reconstructing mostly migrated-in jets. It is
what unfolding is for, but it does mean the highest-`B` bin will be strongly
correction-driven — check it against the closure (`test_mode=0`) before believing it.

The EEC weight reshapes `B` noticeably (weighted gen mean 0.673 vs unweighted 0.711): the
weight `(pt1·pt2)^n` is largest for a balanced pair, so it pulls towards 0.5. A yield run
(`_noeecw`) in `B` is therefore a genuinely different distribution, not a rescale.

### First B results, before any template fit

Produced 2026-09-15 on the **5 × 0.1 binning**, from full QCD/Pythia (10 blocks) and
bJet/Pythia (9 blocks) MC and a full HardProbes data production (50 condor jobs), all under
`OUT_TAG=upartv2_B` so the nominal files were not touched. Numbers below are `both`.

**Cost of a re-production, measured:** 3 min for QCD, 5 min for bJet, ~10 min for the condor
data run — *not* the ~50-60 min the MC section above quotes, because the inputs were warm in
page cache from earlier passes. Iterating on the B binning is cheap when the inputs are hot. Plots and per-bin tables: `plot_z_first_look()` in `template_fit.cpp`, writing to
`results/B_first_look_qcd_pythia_upartv2_B/`.

```bash
root -l -b -q -e '.L template_fit.cpp+' -e 'plot_z_first_look("qcd","pythia","upartv2_B")'
```

**Nothing here is a measurement** — reco level, pre-fit, un-unfolded, no corrections, no SF.

**1. The flavour templates do NOT separate in B.** Unit-area B shapes, 2b / 1b / 0b:

```
0.5-0.6   0.312  0.328  0.301
0.6-0.7   0.286  0.282  0.285
0.7-0.8   0.222  0.212  0.232
0.8-0.9   0.138  0.132  0.135
0.9-1.0   0.042  0.047  0.046
```

All three track each other to within a few percent, and `S/(S+B)` moves only 0.365 → 0.416
across the five bins. **This is good news, not bad:** the background subtraction has almost no
B-dependent shape to impose, so the fit cannot sculpt the B distribution much. It also means
a B-dependent signal fraction buys little over a constant one — worth remembering if the
per-bin fits turn out to be statistics-limited.

**2. m_2B still separates cleanly inside every B bin,** which is what decides whether the
fit works: 2b peaks near 4–5 GeV, 1b/0b near 2–3 GeV, in all five bins. Even the thinnest
bin (0.9–1.0) carries 2565 weighted 2b entries and 28946 effective MC entries.
**The template fit in B will work.**

**3. The response is better in B than in ΔR, and best where the migration is worst.**

```
purity      0.906  0.908  0.913  0.916  0.922
efficiency  0.926  0.931  0.935  0.936  0.942
ΔR purity   0.874 0.908 0.924 0.923 0.923 0.917 0.910 0.894 0.777
```

ΔR's last bin drops to 0.777; B has no weak bin on this measure — it rises monotonically,
and 0.9–1.0 is the *best* bin.

⚠️ **That is not a contradiction of the 0.9–1.0 migration warning above, and the two must not
be conflated.** These are the **analysis's** purity/efficiency — the fraction of reco-passing
jets that also pass the gen selection. The binning study's 0.251 is the **migration diagonal
fraction** — how many jets stay in their own bin. The top B bin is excellent at the first and
poor at the second: almost every jet in it is a genuine selected jet, but many of them
belong in a different B bin. Same word, different quantity.

**4. Data/MC in B rises with B, pre-fit** — monotonic across all five bins:

```
0.939   0.968   1.028   1.109   1.172
```

Now reaching **+17%** in the top bin, a cleaner and stronger trend than the 4-bin version
showed. On the earlier binning it was 0.938→1.105 for `qcd` and 0.934→1.103 for `both`:
**adding the b-enriched sample does not move it**, so the trend is not an artefact of the
QCD sample's flavour composition.

An 18% trend across the range. **It is not the missing UParT SF.** That was worth checking,
because the ΔR SF runs 0.933 → 1.149, a similar size — but ΔR and B turn out to be almost
uncorrelated (EEC-weighted Pearson **+0.049**), and ⟨ΔR⟩ per B bin moves only 0.192 → 0.202,
about a fifth of one ΔR bin. The SF can account for well under 1% of the trend, not 18%.
So the shape difference is real and still to be explained — but it is pre-fit and includes
the 1b/0b background, so the template fit is the next thing that can speak to it, not a
physics conclusion yet.

### ⚠️ B is measured as a YIELD (dN/dB), not EEC-weighted — and why

**The EEC weight is degenerate with the momentum balance.** For a pair with total transverse
momentum `S = pT_b1 + pT_b2`,

```
pT_b1 = B·S,  pT_b2 = (1-B)·S    =>    weight = pT_b1·pT_b2 = B(1-B)·S²
```

so EEC-weighting the B distribution multiplies it by an analytic function **of B itself**.
Measured in data (pT > 100), the mean weight per B bin against that prediction:

```
B bin        <pT_b1·pT_b2>   ratio to bin1   B(1-B) ratio   implied <S²>
0.50-0.60         766           1.000           1.000           3096
0.60-0.70         704           0.919           0.919           3096
0.70-0.80         587           0.767           0.758           3133
0.80-0.90         406           0.529           0.515           3181
0.90-1.00         203           0.265           0.192           4268
```

`<S²>` is constant to ~3% over the first four bins: the weighting **is** `B(1-B) × const`.
It therefore carries essentially no information beyond a kinematic factor, while suppressing
high B by a factor ~5 — exactly the asymmetric-splitting region the observable exists to
probe. (The last bin breaks the pattern, but it is the same bin already flagged for
migration and bias, so it is not evidence of anything.)

**This is the structural difference from ΔR.** There, `pT_b1·pT_b2` is independent of ΔR, so
the energy weight is genuine extra information — energy weighting × angular correlation,
which is what makes an EEC an EEC. For B the weight and the observable are the same
variable, and the correlator degenerates.

**So the B measurement uses the `_noeecw` yield mode of the chain** (`EEC_WEIGHT_OFF=true`),
giving dN/dB. Everything must be unweighted together — templates, data, response, purity,
efficiency, corrections — for the same reason the dR yield run needs its own data
production: EEC-weighted data fitted against unweighted templates mixes two observables.

The EEC-weighted B histograms still exist (`h3D_*_B` in the `_upartv2_B` production) and the
EEC-weighted B fit was run; keep them as the comparison, not as the measurement.

### The dN/dB chain, step by step

```bash
# 1. production (EEC weight OFF), ~7 min for both MC samples with warm inputs
SAMPLE=qcd  GENERATOR=pythia OUT_TAG=upartv2_B EEC_WEIGHT_OFF=true ./run_agg_ntuple_chunks.sh
SAMPLE=bjet GENERATOR=pythia OUT_TAG=upartv2_B EEC_WEIGHT_OFF=true ./run_agg_ntuple_chunks.sh
OUT_TAG=upartv2_B EEC_WEIGHT_OFF=true ./make_hardprobes_condor_scripts.sh && condor_submit condor.submit

# 2. hadd the blocks to the top level -- template_fit.cpp reads MERGED files, not per-block
hadd -f <...>_fMCGEN_noeecw_upartv2_B.root <...>/block_*/<same name>.root

# 3. the fit
root -l -b -q 'template_fit.cpp+("both","pythia",false,true,"B")'
#                                                      ^^^^ EEC_WEIGHT_OFF
```

Sanity check on step 1/2: with the weight off, MC `h3D_*_B` must equal `h_count_*_B`
**exactly** — verified, 0 differing cells.

⚠️ **Data is the exception, and this is not a bug.** `h3D_data_B` and `h_count_data_B`
differ (1.949M vs 1.773M) because the **trigger prescale multiplies `eec`, never the
count** — so with the weight off the fill weight is `1 × prescale`. The yield measurement
wants the prescale-corrected `h3D_data_B`; `h_count_data_B` is the raw number of jets
recorded. The fit reads `h3D_data*`, so it gets the right one.

**The fitted dN/dB signal fractions**, against the EEC-weighted ones:

```
                  B bin      yield    err      EEC-wtd    diff
80 < pT < 100   0.5-0.6      0.3903  0.0027    0.3714    +0.019
                0.6-0.7      0.3811  0.0027    0.3605    +0.021
                0.7-0.8      0.3754  0.0028    0.3564    +0.019
                0.8-0.9      0.3942  0.0033    0.3707    +0.024
                0.9-1.0      0.3717  0.0047    0.3303    +0.041
pT > 100        0.5-0.6      0.4847  0.0023    0.4712    +0.014
                0.6-0.7      0.4851  0.0023    0.4769    +0.008
                0.7-0.8      0.4810  0.0023    0.4784    +0.003
                0.8-0.9      0.4787  0.0026    0.4787    -0.000
                0.9-1.0      0.4678  0.0035    0.4662    +0.002
```

Still flat in B, as expected — the flavour templates do not separate in B either way. The
yield fractions sit slightly **above** the EEC-weighted ones and are slightly more precise
(the EEC weight adds variance without adding information here). The gap is largest in the
lower pT bin and at high B, i.e. exactly where the `B(1-B)` weighting bites hardest.

Output: `TemplateFits_both_pythia_noeecw_B_upartv2/`.

### Where the raw yields are

The `h_count_*` histograms — **not** the `h3D_*` ones. Same jets, same selection, same
binning; the only difference is that `h_count_*` is filled with the tree weight alone and
`h3D_*` carries the EEC weight on top. That is the difference between "how many jets" and
"how much EEC", and it is why an EEC-weighted plot is never a yield.

| | data file | MC file |
|---|---|---|
| ΔR | `h_count_data` | `h_count_bb`, `h_count_b`, `h_count_0b` |
| B | `h_count_data_B` | `h_count_bb_B`, `h_count_b_B`, `h_count_0b_B` |

Section 5 of `plot_z_first_look()` prints and plots them.

**Data, raw b-tagged 2-SV jet counts per B bin** (full HardProbes, 5 × 0.1 binning):

```
0.5-0.6     389774    sqrt(N) 624    0.16%
0.6-0.7     385502            621    0.16%
0.7-0.8     377748            615    0.16%
0.8-0.9     367653            606    0.16%
0.9-1.0     252410            502    0.20%
total      1773087 jets
```

⚠️ The trigger **prescale is not applied** in `h_count_*` — the prescale multiplies `eec`,
not the count — so these are jets actually recorded, not a luminosity-scaled yield.

⚠️ The counts are **not** flat per bin even after dividing by bin width (3.89, 3.88, 3.81,
3.21 ×10⁶ per unit B). That is expected: the binning was chosen equal-occupancy on the
**EEC-weighted gen** distribution, not on unweighted reco data counts.

**MC template statistics — the number that matters for the fit.** MC `h_count_*` are
*weighted* predictions, so the useful quantity is effective entries, `(Σw)²/Σw²`, per bin:

```
2b effective entries, both:   55614   50352   50690   44661   28946    (~0.4-0.6% per bin)
```

Running bjet (2026-09-15) bought a **factor ~15** in template statistics, exactly as the
"bjet is there for template statistics" note predicts. With `qcd` alone the per-B-bin fit was
limited by MC templates (~2%) against data at 0.16–0.20%; with `both` the templates are
~0.5% and the gap is much narrower. **Fit the B templates on `both`.**

⚠️ But `S/(S+B)` and the **0b** column move when bjet is added (0b weighted counts fall from
~0.15–0.31 to ~0.26–0.55 against a 2b column that grew 15×, i.e. the mistag *fraction*
collapses), because the b-enriched sample contains almost no mistagged light jets. This is
the same trap `template_fit.cpp` already guards against: `c' = 0B/(0B+1B)` is computed from
the **qcd integrals only**, never the bjet ones, precisely because bjet would bias it
downward. When the B fit is written, `c'` must keep that rule — take the shapes from `both`,
take `c'` from `qcd`.

### The B template fit

```bash
root -l -b -q 'template_fit.cpp+("both","pythia",false,false,"B")'
```

The 5th argument is the observable, `dr` (default) or `B`. `dr` is unchanged in every
respect — same inputs, same histogram names, same output directory — **verified**: all four
variations (`nominal`, `var0B_2`, `var0B_0`, `fitRange0to8`) come out bit-identical, 176
histograms each, max |Δ| = 0.

How the axis was generalised: `gFitObs()` in `observables.h` is a driver-set global, exactly
like `sDirname`, read by `CheckInputBinning()`, the histogram names and the plot labels. It
defaults to dR, so any path that does not set it behaves as before.

⚠️ **The B fit runs on a DIFFERENT production from the nominal dR fit.** The nominal dR fit
reads Afnan's pre-merged Pythia files (`/data_CMS/cms/shatat/…MattProd/…merged.root`) and
her data file; those predate the observable axis and contain **no B histograms**. The B fit
therefore reads the local `_upartv2_B` production — MC hadd'ed to the `agg_ntuple_chunks` top
level, data to `agg_template_chunks`. That set is internally consistent (its own templates
and its own data), but a B signal fraction and the nominal dR one are **not** from the same
MC. To compare them like for like, the dR fit would have to be re-run on `_upartv2_B` too —
`mcgenTemplates()` does not offer that yet.

Every B fit reads a `_upartv2_B` production, whatever the generator or variation:
`mcgenTemplates()` builds the path from the same tags `apply_unfolding_2d.C` reads its blocks
with (`<sample>_fMCGEN[_trkdrop030][_noeecw]_upartv2_B.root`, under `QCD`/`bJet` or
`QCDHerwig`/`bJetHerwig`). It never falls back to a ΔR file. A combination that was never
produced stops on the missing-file check, naming the file. (Until 2026-09-24 Herwig and the
track drop were refused outright for B; that guard is kept commented in the driver.)

Output: `TemplateFits_both_pythia_B_upartv2/` — the `_B` tag keeps it out of the dR fit's
folder, which matters because `apply_unfolding_2d.C` reads `h_sig_fraction_fit` **by name**
and a B-binned one sitting there would be silently wrong.

**The result: the B signal fraction is flat.**

```
             0.5-0.6  0.6-0.7  0.7-0.8  0.8-0.9  0.9-1.0
80-100 GeV    0.371    0.361    0.356    0.371    0.330
100-120 GeV   0.471    0.477    0.478    0.479    0.466
errors        ~0.003                              ~0.005
```

against the dR fit, which swings **0.124 → 0.930** across its nine bins.

This is the first-look prediction confirmed at the fit level: the flavour templates do not
separate in B, so the background subtraction in B is close to a **constant rescale rather
than a shape distortion**. Practically, the B measurement is much less exposed to the fit
than the dR one is — the mistag and template-modeling systematics, which enter through the
signal fraction's *shape*, have very little shape to move here.

The last B bin (0.9–1.0) has both the lowest fraction and the largest error, as expected
from its statistics; it is the same bin flagged above for migration.

### The scale-factor decision

**The `B` measurement runs without the UParT scale factor, by choice (2026-09-15).**
`lifetime_jp_sfb_agg_dr_wp712_eec_systematics.root` is binned in reco **ΔR** over the 9
analysis bins, and there is no `B`-binned equivalent, so there is nothing to apply. For ΔR
the SF is a *correction, not an option* and every data run applies it; for `B` it is simply
absent. Read a `B` result accordingly — it is **not** the ΔR chain minus a systematic, it is
missing a correction the ΔR result has, and the ΔR SF runs 0.933 to 1.149 across its range,
so the size of what is missing is not negligible. Getting it later means either applying the
SF in 2D at fill time (both ΔR and `B` are available on the same jet, so this is possible
without a new calibration) or a `B`-binned measurement of it.

### The full B band — every input, in order

The code is ready for all of it (2026-09-24); what it needs is **six MC productions, one data
production, seven hadds, three fits and seven unfoldings**, all under `OUT_TAG=upartv2_B` with
the EEC weight off. Everything below runs on LLR in `workflow/` after
`source setup_roounfold_env.sh`.

> ⚠️ **Blocked on the input chunks (found 2026-09-24).** Marc's `…_UParTV2_chunks` areas were
> re-made between 2026-09-22 and 2026-09-24 as `merged_block_NNNN_…_UParTV2_rho.root`. The
> old MC inputs are **gone**, and there are now 8 Pythia blocks per sample where there were
> 10 (QCD) and 9 (bJet). The data re-merge was **still being written** on the morning of
> 2026-09-24 (`HardProbes/<N>/UParTV2_chunks/merged_block_000[3-9]_UParTV2.root` were 543-byte
> stubs). Two consequences:
> - `run_agg_ntuple_chunks.sh` globs `merged_block_*_<FILE_TAG>.root` and finds **nothing**;
>   `make_hardprobes_condor_scripts.sh` points at the **stubs**, which pass its
>   missing-input check. Neither is safe to run until the inputs settle.
> - Any B production made now comes from **different MC than the EEC(ΔR) result**, which
>   was made from the old chunks. Decide whether ΔR is re-produced on the `_rho` inputs too
>   before spending the runs.

```bash
# 1. MC -- ONE AT A TIME: each run rm -rf's its staging dirs and recompiles the shared .so.
#    Nominal Pythia (the measurement and the 0B pair), Herwig (template modelling AND detector
#    response), then the 3% track drop.
for g in pythia herwig; do for s in qcd bjet; do
  SAMPLE=$s GENERATOR=$g OUT_TAG=upartv2_B EEC_WEIGHT_OFF=true ./run_agg_ntuple_chunks.sh
done; done
for s in qcd bjet; do
  SAMPLE=$s GENERATOR=pythia OUT_TAG=upartv2_B EEC_WEIGHT_OFF=true TRACK_EFF_UNC=true ./run_agg_ntuple_chunks.sh
done

# 2. Data -- AFTER step 1: condor jobs load the .so each MC run deletes and rebuilds.
#    Data is never varied, so one production serves every entry of the band.
OUT_TAG=upartv2_B EEC_WEIGHT_OFF=true ./make_hardprobes_condor_scripts.sh && condor_submit condor.submit

# 3. hadd -- template_fit.cpp reads MERGED MCGEN files; apply_unfolding_2d.C reads the blocks.
P=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024
for d in QCD:qcd:_noeecw bJet:bjet:_noeecw QCDHerwig:qcd:_noeecw bJetHerwig:bjet:_noeecw \
         QCD:qcd:_trkdrop030_noeecw bJet:bjet:_trkdrop030_noeecw; do
  IFS=: read sub s v <<< "$d"
  f=Run3_btagWP0712_template_for_fit_histos_3D_${s}_fMCGEN${v}_upartv2_B.root
  hadd -f $P/$sub/agg_ntuple_chunks/$f $P/$sub/agg_ntuple_chunks/block_*/$f
done
f=Run3_btagWP0712_template_for_fit_histos_3D_data_fMCGEN_noeecw_upartv2_B.root
hadd -f $P/HardProbes/agg_template_chunks/$f $P/HardProbes/agg_template_chunks/HardProbes*/block_*/$f

# 4. Fits. The nominal one also writes the var0B_2 / var0B_0 refits the mistag pair reads.
root -l -b -q 'template_fit.cpp+("both","pythia",false,true,"B")'
root -l -b -q 'template_fit.cpp+("both","pythia",true, true,"B")'   # tracking_eff
root -l -b -q 'template_fit.cpp+("both","herwig",false,true,"B")'   # template_fit

# 5. Unfoldings. Arguments: SAMPLE, UNFOLDING_GEN, test_mode, bayes, scan, TF_GEN, TRACK_EFF_UNC,
#    TF_VARIATION, SFUPART ("off": none exists for B), EEC_WEIGHT_OFF, OBSERVABLE.
u() { root -l -b -q "apply_unfolding_2d.C(\"both\",\"$1\",$2,false,false,\"$3\",$4,\"$5\",\"off\",true,\"B\")"; }
u pythia 0 pythia false nominal   # full-MC closure -- must be 1.0000 in every bin
u herwig 0 pythia false nominal   # the same for the Herwig MC
u pythia 2 pythia false nominal   # THE RESULT
u pythia 2 herwig false nominal   # template_fit     (MC template modeling)
u herwig 2 pythia false nominal   # unfolding_model  (Detector response)
u pythia 2 pythia true  nominal   # tracking_eff     (Tracking efficiency)
u pythia 2 pythia false var0B_2   # mistag_0B_up     (Light jet mistagging)
u pythia 2 pythia false var0B_0   # mistag_0B_down

# 6. The band, data_vs_gen and systematics_curves, into the nominal B folder.
root -l -b -q 'apply_weights_and_systematics.C("both","pythia",2,false,"h_data_fully_corrected_2D","pythia",2,true,"B")'
```

Timing on warm inputs was ~3 min (QCD) and ~5 min (bJet) per Pythia run on 2026-09-15. Herwig
has never been run with the B axis and the `_rho` inputs are new, so expect the cold-cache
~50–60 min per sample instead.

## Conventions

- Compile the `.so` **once**; jobs `gSystem->Load()` it. Parallel ACLiC compilation in the
  same directory produces half-written artifacts.
- ACLiC products, ROOT files, plots and logs are gitignored — do not commit them.
- `*~`, `#*#` (Emacs) and `._*` (macOS/sshfs) files litter the tree; ignore them.
- Unwanted analysis code is commented out with `/* ---- disabled (kept for reference) ---- */`
  rather than deleted.
