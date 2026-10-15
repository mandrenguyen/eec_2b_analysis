#!/bin/bash

set -u

WORK=/home/llr/cms/zaidan/analysis_lise/eec_2b_analysis/workflow
OUT_BASE=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/HardProbes/agg_template_chunks

# ACLIC_BUILD_DIR must match the build dir set in rootlogon.C.
ACLIC_BUILD_DIR=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/build
SCRIPT_DIR=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/condor/hardprobes_data_scripts
CONDOR_LOG_DIR=/data_CMS/cms/zaidan/bJetAggRun3/PPRef2024/condor/logfiles

# INPUT_TAG names the chunk subdir (<TAG>_chunks) and the per-block filenames.
#INPUT_TAG=negTagFix
INPUT_TAG=UParTV2
# Marc re-made the data chunks on 2026-09-22..24 as merged_block_NNNN_UParTV2_rho.root and the
# un-suffixed ones are gone (they were 543-byte stubs mid-rewrite -- which the missing-input
# check below lets through, so never point this at a half-written area).
INPUT_SUFFIX=${INPUT_SUFFIX-_rho}

# SAMPLE_TAG replaces the old dataType argument: the macro puts it in the output histogram
# names (..._template_for_fit_histos_3D_<SAMPLE_TAG>_f). These jobs run HardProbes data.
SAMPLE_TAG=data

# ---------------------------------------------------------------------------
# EEC weight on/off -- the OBSERVABLE, not a systematic.
# ---------------------------------------------------------------------------
# false (default) measures the EEC: every 2b jet enters weighted by (pt1*pt2)^n.
# true drops that weight, so h3D_data becomes the plain 2b-jet yield per dr bin.
# DATA MUST MATCH THE MC: unfolding EEC-weighted data through an unweighted
# response (or fitting it against unweighted templates) mixes two observables
# and is simply wrong, so a yield analysis needs this production as well as the
# MC one. Overridable for one run:  EEC_WEIGHT_OFF=true ./make_hardprobes_condor_scripts.sh
EEC_WEIGHT_OFF=${EEC_WEIGHT_OFF:-false}   # true | false
EECW_TAG=$([ "${EEC_WEIGHT_OFF}" = true ] && echo "_noeecw" || echo "")
if [ "${EEC_WEIGHT_OFF}" != true ] && [ "${EEC_WEIGHT_OFF}" != false ]; then
  echo "EEC_WEIGHT_OFF must be exactly 'true' or 'false', got '${EEC_WEIGHT_OFF}'"; exit 1
fi

# Suffix appended on the way out of the staging dir. The MACRO already writes the
# observable tag into its own filename (EecWeight::tag()), exactly as it writes the
# tracking tag, so EECW_TAG must NOT be repeated here -- doing so produced
# ..._fMCGEN_noeecw_noeecw_upartv2.root. Same rule as run_agg_ntuple_chunks.sh: the
# variation tag names the staging dir and the logs, the macro names the file.
# Overridable, same as in run_agg_ntuple_chunks.sh, so an exploratory data production can
# write to its own files:  OUT_TAG=upartv2_B ./make_hardprobes_condor_scripts.sh
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
# Staging dir only, so two productions can never share a staging area.
STAGE_TAG="${EECW_TAG}${OUT_SUFFIX}"
# The job scripts get their own directory per production too, so two productions (the EEC and
# the yield one) can be generated and queued at the same time without the second rm'ing the
# first one's scripts. Submit with the matching Executable (printed at the end).
SCRIPT_DIR="${SCRIPT_DIR}${STAGE_TAG}"

# Must match BTAG_WP in run_agg_ntuple_chunks.sh.
case "${INPUT_TAG}" in
  UParTV2)    BTAG_WP=0.712 ;;
  negTagFix)  BTAG_WP=0.868 ;;
  *) echo "unknown INPUT_TAG '${INPUT_TAG}': set BTAG_WP for it explicitly"; exit 1 ;;
esac
echo "sample ${INPUT_TAG}, b-tag WP ${BTAG_WP}"
if [ "${EEC_WEIGHT_OFF}" = true ]; then
  echo "EEC WEIGHT OFF: this data production measures YIELDS (dN/dr), file tag ${EECW_TAG}"
fi

mkdir -p "${SCRIPT_DIR}" "${CONDOR_LOG_DIR}"
rm -f "${SCRIPT_DIR}"/job_*.sh "${SCRIPT_DIR}"/jobs_*.sh

job_index=0
for primary_dataset in $(seq 0 4); do
  input_dir=/data_CMS/cms/mnguyen/bJetAggRun3/PPRef2024/HardProbes/${primary_dataset}/${INPUT_TAG}_chunks

  for i in $(seq 0 9); do
    block=$(printf "000%d" "${i}")
    script="${SCRIPT_DIR}/jobs_${job_index}.sh"
    input="${input_dir}/merged_block_${block}_${INPUT_TAG}${INPUT_SUFFIX}.root"
    outdir="${OUT_BASE}/HardProbes${primary_dataset}/block_${block}"
    logdir="${OUT_BASE}/HardProbes${primary_dataset}/logs"
    # The macro picks its own filenames, so stage them and rename on the way out.
    stagedir="${outdir}/.stage${STAGE_TAG}"

    cat > "${script}" <<EOF
#!/bin/bash

set -euo pipefail

WORK=${WORK}
JOB_INDEX=${job_index}
PRIMARY_DATASET=${primary_dataset}
BLOCK=${block}
INPUT=${input}
OUTDIR=${outdir}
LOGDIR=${logdir}
STAGEDIR=${stagedir}
OUT_SUFFIX=${OUT_SUFFIX}
ACLIC_BUILD_DIR=${ACLIC_BUILD_DIR}

cd "\${WORK}"
source "\${WORK}/setup_roounfold_env.sh"
mkdir -p "\${OUTDIR}" "\${LOGDIR}"
rm -rf "\${STAGEDIR}"
mkdir -p "\${STAGEDIR}"

if [ ! -f "\${INPUT}" ]; then
  echo "ERROR: missing input \${INPUT}"
  exit 2
fi

if [ ! -f "\${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so" ]; then
  echo "ERROR: missing \${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so -- run ./run_agg_ntuple_chunks.sh style compile (LCG env) once before submitting"
  exit 3
fi

root -l -b -q -e "gSystem->AddIncludePath(\"-I\${ROOUNFOLD_INC} -I\${ROOUNFOLD_BUILD}\"); if (gSystem->Load(\"\${ROOUNFOLD_BUILD}/libRooUnfold.so\") < 0) gSystem->Exit(3); if (gSystem->Load(\"\${ACLIC_BUILD_DIR}/create_files_for_template_fit_cpp.so\") < 0) gSystem->Exit(3); create_files_for_template_fit(3,80,2,1,true,false,${BTAG_WP},true,false,true,0,-1,\"\${INPUT}\",\"\${STAGEDIR}\",\"${SAMPLE_TAG}\",false,20260908,${EEC_WEIGHT_OFF},${APPLY_JEC})"

for f in "\${STAGEDIR}"/*.root; do
  [ -e "\${f}" ] || continue
  base=\$(basename "\${f}")
  mv "\${f}" "\${OUTDIR}/\${base%.root}\${OUT_SUFFIX}.root"
done
rmdir "\${STAGEDIR}" 2>/dev/null || true
EOF

    chmod +x "${script}"
    job_index=$((job_index + 1))
  done
done

echo "created 50 scripts in ${SCRIPT_DIR}"
echo "submit with:"
echo "  condor_submit condor.submit -append 'Executable = ${SCRIPT_DIR}/jobs_\$(Process).sh'"
