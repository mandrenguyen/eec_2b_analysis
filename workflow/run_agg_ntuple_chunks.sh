#!/bin/bash
set -u


SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK=${SCRIPT_DIR}
source "${WORK}/setup_roounfold_env.sh"


# ============================ SET THE SAMPLE HERE ============================
# Two flags; input dir, output dir, chunk filenames and the macro's sample tag all
# follow from them.
#
# Each knob below is written VAR=${VAR:-default}, so the value here is the default and
# editing it still works exactly as before -- but a run can also override it for one
# invocation without touching the file, which is what makes a two-sample production a
# loop instead of two edits:
#     SAMPLE=qcd EEC_WEIGHT_OFF=true ./run_agg_ntuple_chunks.sh
SAMPLE=${SAMPLE:-bjet}          # qcd | bjet
GENERATOR=${GENERATOR:-pythia}  # pythia | herwig

# b-tag production of the chunks. Data runs over negTagFix, so use negTagFix here too;
# negTag is the older production. Only Pythia8 QCD has negTag/negTagFix chunks.
INPUT_TAG=UParTV2   # UParTV2 | negTagFix | negTag

case "${SAMPLE}_${GENERATOR}" in
  qcd_pythia)   IN_SUBDIR="QCD/Pythia8_${INPUT_TAG}_chunks";  OUT_SUBDIR=QCD        ;;
  bjet_pythia)  IN_SUBDIR="bJet/Pythia8_${INPUT_TAG}_chunks"; OUT_SUBDIR=bJet       ;;
  qcd_herwig)   IN_SUBDIR="QCDHerwig/${INPUT_TAG}_chunks";    OUT_SUBDIR=QCDHerwig  ;;
  bjet_herwig)  IN_SUBDIR="bJetHerwig/${INPUT_TAG}_chunks";   OUT_SUBDIR=bJetHerwig ;;
  *) echo "unknown SAMPLE/GENERATOR '${SAMPLE}/${GENERATOR}': use qcd|bjet and pythia|herwig"; exit 1 ;;
esac

# Pythia8 chunks carry the generator in the filename, the Herwig ones do not:
#   merged_block_0000_Pythia8_UParTV2.root  vs  merged_block_0000_UParTV2.root
# INPUT_SUFFIX: Marc re-made every chunk on 2026-09-22..24 as ..._UParTV2_rho.root and the
# un-suffixed ones are gone. Keep OUT_TAG distinct from any pre-_rho production (upartv2rho,
# upartv2rho_B -- see eraTag() in result_paths.h): _rho has 8 Pythia blocks per sample where
# the old had 10/9, and old blocks left under the same name would be summed in.
INPUT_SUFFIX=${INPUT_SUFFIX-_rho}
FILE_TAG=$([ "${GENERATOR}" = pythia ] && echo "Pythia8_${INPUT_TAG}" || echo "${INPUT_TAG}")${INPUT_SUFFIX}

INPUT_DIR=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/${IN_SUBDIR}
OUT_BASE=$mydata/bJetAggRun3/PPRef2024/${OUT_SUBDIR}/agg_ntuple_chunks


SAMPLE_TAG=${SAMPLE}

# Overridable like the knobs above, so an exploratory production can write to its own files
# instead of overwriting the nominal one:  OUT_TAG=upartv2_B ./run_agg_ntuple_chunks.sh
OUT_TAG=${OUT_TAG:-upartv2}
OUT_SUFFIX="${OUT_TAG:+_${OUT_TAG}}"

# ---------------------------------------------------------------------------
# Jet energy corrections (2024 pp-reference JEC, as in Matt's JEC_Prompt24HIpp files)
# ---------------------------------------------------------------------------
# false (default) keeps jtpt as stored in the forest, so every existing production is
# reproduced unchanged. true recomputes jtpt from rawpt with L1FastJet + L2Relative (+ the
# L2Residual on data) before any selection -- see the Jec namespace in
# create_files_for_template_fit.cpp. The OUT_TAG must then say "jec" (eraTag() in
# result_paths.h: upartv2rhojec / upartv2rhojec_B), and must NOT say it otherwise, so
# corrected and uncorrected files can never share a name.
APPLY_JEC=${APPLY_JEC:-false}   # true | false
if [ "${APPLY_JEC}" != true ] && [ "${APPLY_JEC}" != false ]; then
  echo "APPLY_JEC must be exactly 'true' or 'false', got '${APPLY_JEC}'"; exit 1
fi
case "${OUT_TAG}" in
  *jec*) [ "${APPLY_JEC}" = true ]  || { echo "OUT_TAG '${OUT_TAG}' says jec but APPLY_JEC=false -- refusing"; exit 1; } ;;
  *)     [ "${APPLY_JEC}" = false ] || { echo "APPLY_JEC=true needs an OUT_TAG containing 'jec' (e.g. upartv2rhojec) -- refusing"; exit 1; } ;;
esac
[ "${APPLY_JEC}" = true ] && echo "JEC ON: jtpt recomputed from rawpt (L1FastJet + L2Relative, + L2Residual on data)"

# ---- Tracking-efficiency systematic -----------------------------------------
# false = nominal production (nothing is dropped, output names unchanged)
# true  = the systematic: throw away 3% of the RECONSTRUCTED tracks in MC during
#         the B reconstruction, then redo the unfolding with these files and take
#         the symmetrised difference against the nominal. The 3% lives in
#         TrkEffSyst::kDropFraction, not here.
#
# Which tracks go is a hash of (TRK_SEED, entry number, track index), not a draw
# from a stream, so this is reproducible: rerun the same block, rerun with a
# different number of blocks, or rerun serially instead of in parallel and the
# same tracks are dropped every time. Change TRK_SEED only if you deliberately
# want an independent throw. See tracking_efficiency_syst.h.
#
# The macro tags its own output files (_trkdrop030 for 3%), so a variation run can
# never overwrite the nominal one and both can sit in the same block directory.
TRACK_EFF_UNC=${TRACK_EFF_UNC:-false}   # true | false

TRK_SEED=20260908

# ---------------------------------------------------------------------------
# EEC weight on/off -- this is the OBSERVABLE, not a systematic.
# ---------------------------------------------------------------------------
# false (default) measures the EEC: every 2b jet enters weighted by (pt1*pt2)^n.
# true drops that weight everywhere -- templates, response, purity, efficiency,
# truth -- so the same chain measures YIELDS and the unfolded result is dN/dr.
# Output is tagged "_noeecw", so a yield production can never overwrite an EEC
# one and both can live in the same block directory. See the EecWeight namespace
# at the top of create_files_for_template_fit.cpp.
EEC_WEIGHT_OFF=${EEC_WEIGHT_OFF:-false}   # true | false

# ---------------------------------------------------------------------------
# SV b scale factor (sv_sfb_singleb_aggpt_jp_wp712_mu5_systematics.root, see sv_sf.h)
# ---------------------------------------------------------------------------
# off     no SF, output names as before 2026-10-12 (was the default until then)
# nominal every b-matched reco SV weighted by SF(aggregate pT)           tag _svsf
# up/down SF shifted by +/- sqrt(stat^2 + syst_total^2), coherently    tag _svsfup / _svsfdn
# MC only. The macro tags its own files, as for the tracking variation.
# Default "nominal" since 2026-10-12, matching the readers (svsfNominal() in result_paths.h):
# the SF is part of the nominal MC now. SV_SF=off reproduces the earlier no-SF productions.
# SV_SF=${SV_SF:-off}   # off | nominal | up | down
SV_SF=${SV_SF:-nominal}   # off | nominal | up | down
case "${SV_SF}" in
  off)     SVSF_TAG="" ;;
  nominal) SVSF_TAG="_svsf" ;;
  up)      SVSF_TAG="_svsfup" ;;
  down)    SVSF_TAG="_svsfdn" ;;
  # per-contribution variations (2026-10-13), see sv_sf.h
  statup|statdn|nomuonup|nomuondn|bbup|bbdn|jpup|jpdn|qqup|qqdn) SVSF_TAG="_svsf${SV_SF}" ;;
  *) echo "SV_SF must be off | nominal | up | down | stat|nomuon|bb|jp|qq + up|dn, got '${SV_SF}'"; exit 1 ;;
esac
[ "${SV_SF}" != off ] && echo "SV b SCALE FACTOR: ${SV_SF}, file tag ${SVSF_TAG}"

# ---------------------------------------------------------------------------
# Track-BDT efficiency and mistag scale factors
# ---------------------------------------------------------------------------
# nominal applies SFb and SFl through the four category-dependent effective WPs.
# boff/loff remove one coupled in/out correction at a time; off is the uncorrected
# cross-check. This is MC-only; the data production keeps the 0.365 working point.
TRACK_BDT_SF=${TRACK_BDT_SF:-nominal}   # off | nominal | boff | loff
case "${TRACK_BDT_SF}" in
  off)     TRACKBDT_TAG="" ;;
  nominal) TRACKBDT_TAG="_trkbdtsf" ;;
  boff)    TRACKBDT_TAG="_trkbdtsf_boff" ;;
  loff)    TRACKBDT_TAG="_trkbdtsf_loff" ;;
  *) echo "TRACK_BDT_SF must be off | nominal | boff | loff, got '${TRACK_BDT_SF}'"; exit 1 ;;
esac
echo "Track-BDT SF: ${TRACK_BDT_SF}${TRACKBDT_TAG:+, file tag ${TRACKBDT_TAG}}"

# Same tags the macro builds, for the staging dir and the log filenames.
TRK_TAG=$([ "${TRACK_EFF_UNC}" = true ] && echo "_trkdrop030" || echo "")
EECW_TAG=$([ "${EEC_WEIGHT_OFF}" = true ] && echo "_noeecw" || echo "")
TRK_TAG="${TRK_TAG}${SVSF_TAG}${TRACKBDT_TAG}${EECW_TAG}"

case "${INPUT_TAG}" in
  UParTV2)   BTAG_WP=0.712 ;;
  negTagFix) BTAG_WP=0.868 ;;
  negTag)    BTAG_WP=0.868 ;;
  *) echo "unknown INPUT_TAG '${INPUT_TAG}': set BTAG_WP for it explicitly"; exit 1 ;;
esac

echo "sample ${SAMPLE} ${GENERATOR} (${INPUT_TAG}), b-tag WP ${BTAG_WP}, sample tag ${SAMPLE_TAG}"
if [ "${TRACK_EFF_UNC}" = true ]; then
  echo "TRACKING-EFFICIENCY SYSTEMATIC ON: 3% of reco tracks dropped (seed ${TRK_SEED}), file tag ${TRK_TAG}"
elif [ "${TRACK_EFF_UNC}" = false ]; then
  echo "tracking-efficiency systematic off (nominal)"
else
  echo "TRACK_EFF_UNC must be exactly 'true' or 'false', got '${TRACK_EFF_UNC}'"; exit 1
fi
if [ "${EEC_WEIGHT_OFF}" = true ]; then
  echo "EEC WEIGHT OFF: this run measures YIELDS (dN/dr), not the EEC, file tag ${EECW_TAG}"
elif [ "${EEC_WEIGHT_OFF}" != false ]; then
  echo "EEC_WEIGHT_OFF must be exactly 'true' or 'false', got '${EEC_WEIGHT_OFF}'"; exit 1
fi
echo "input  ${INPUT_DIR}"
echo "output ${OUT_BASE}"
# =============================================================================


LOG_DIR=${OUT_BASE}/logs

ACLIC_BUILD_DIR=$mydata/bJetAggRun3/PPRef2024/build

mkdir -p "${LOG_DIR}"

cd "${WORK}" || exit 1

COMPILE_LOG="${LOG_DIR}/compile.log"
SO="${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so"

# SKIP_COMPILE=true reuses the .so already in the build dir instead of deleting and
# rebuilding it. That is what lets several productions (and the condor data jobs) run at
# the same time: a normal run rm's the .so, and every other job that loads it meanwhile
# fails. The staleness check below still runs, so a .so older than the sources is refused.
SKIP_COMPILE=${SKIP_COMPILE:-false}
if [ "${SKIP_COMPILE}" = true ]; then
  echo "SKIP_COMPILE: reusing ${SO}"
  [ -f "${SO}" ] || { echo "no ${SO} -- run once without SKIP_COMPILE first"; exit 1; }
else
echo "compiling create_files_for_template_fit.cpp once before launching parallel jobs"

# Drop half-written ACLiC products from a previous parallel launch.
rm -f "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp.d \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp.so \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp_ACLiC_dict.* \
      "${ACLIC_BUILD_DIR}"/create_files_for_template_fit_cpp_ACLiC_map.*

## compilation of create_files_for_template_fit.cpp
#
# ⚠️ `root -l -b` EXITS 0 EVEN WHEN ACLiC FAILS. Its status says only that ROOT started and
# quit, not that the library got built. Trusting it launched a full 19-block production
# against a STALE .so on 2026-09-21 -- every job silently measured the PREVIOUS code's
# observables and wrote them under the new OUT_TAG. So the .so itself is what is checked
# below, not the exit code:
#   - the file must exist (the rm above deleted any earlier one, so it can only be there if
#     this compile produced it), and
#   - the compile output must not contain an ACLiC error.
root -l -b > "${COMPILE_LOG}" 2>&1 <<EOF
.L create_files_for_template_fit.cpp++
.q
EOF
compile_status=$?

if [ ${compile_status} -ne 0 ] || [ ! -f "${SO}" ] || grep -q "Error in <ACLiC>" "${COMPILE_LOG}"; then
  echo "compile FAILED -- not launching any jobs. Last 80 lines of ${COMPILE_LOG}:"
  tail -n 80 "${COMPILE_LOG}"
  exit 1
fi
fi   # SKIP_COMPILE

# The .so must be NEWER than every source it is built from, or a failed rebuild left the
# previous one in place and the jobs would run the wrong code.
for src in create_files_for_template_fit.cpp observables.h binning_histos_small.h \
           result_paths.h tTree.h tracking_efficiency_syst.h central_selections.h sv_sf.h \
           track_bdt_sf.h; do
  if [ -e "${WORK}/${src}" ] && [ "${WORK}/${src}" -nt "${SO}" ]; then
    echo "compile FAILED -- ${src} is newer than ${SO}; the build did not take. Not launching."
    tail -n 40 "${COMPILE_LOG}"
    exit 1
  fi
done

echo "compile finished; launching chunk jobs"



shopt -s nullglob
inputs=("${INPUT_DIR}"/merged_block_*_${FILE_TAG}.root)
shopt -u nullglob

if [ ${#inputs[@]} -eq 0 ]; then
  echo "no chunks matching ${INPUT_DIR}/merged_block_*_${FILE_TAG}.root"
  exit 1
fi
echo "found ${#inputs[@]} chunks"


for input in "${inputs[@]}"; do

  base=$(basename "${input}")
  block=${base#merged_block_}
  block=${block%%_*}

  outdir="${OUT_BASE}/block_${block}"
  stagedir="${outdir}/.stage${TRK_TAG}${OUT_SUFFIX}"
  rm -rf "${stagedir}"
  mkdir -p "${stagedir}"

  (
  
  # Remeber: create_files_for_template_fit(..., bool apply_jec = false, const char* sv_sf = "off", const char* track_bdt_sf = "nominal")
  # MC vs data is the isMC argument; bjet vs qcd is SAMPLE_TAG.

  nice -n 10 root -l -b -q -e "gSystem->Load(\"${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\"); create_files_for_template_fit(3,80,2,1,true,true,${BTAG_WP},true,true,true,0,-1,\"${input}\",\"${stagedir}\",\"${SAMPLE_TAG}\",${TRACK_EFF_UNC},${TRK_SEED},${EEC_WEIGHT_OFF},${APPLY_JEC},\"${SV_SF}\",\"${TRACK_BDT_SF}\")"

    for f in "${stagedir}"/*.root; do
      [ -e "${f}" ] || continue
      base=$(basename "${f}")
      mv "${f}" "${outdir}/${base%.root}${OUT_SUFFIX}.root"
    done
    rmdir "${stagedir}" 2>/dev/null
  ) > "${LOG_DIR}/block_${block}${TRK_TAG}${OUT_SUFFIX}.log" 2>&1 &


  echo "submitted block ${block}, pid $!"
done

wait
echo "all chunk jobs finished"
