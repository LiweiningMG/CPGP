#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# CPGP main workflow
#
# Edit PROJECT_ROOT before running. All other paths are derived from it.
###############################################################################

PROJECT_ROOT="/path/to/CPGP"

DATASETS=("Xie2021" "Lee2019")
THREADS=25
FOLD=5
REP=10
MCMC_CYCLES=30000
BURNIN=20000
THIN=10

# Set RUN_STEPS to choose which stages to run.
# Available stages: within joint bayes mgrm
RUN_STEPS=("within" "joint" "bayes" "mgrm")

if [[ "${PROJECT_ROOT}" == "/path/to/CPGP" ]]; then
  echo "Please edit PROJECT_ROOT in main.sh before running."
  exit 1
fi

CODE_DIR="${PROJECT_ROOT}/code"
BIN_DIR="${CODE_DIR}/bin"
SHELL_DIR="${CODE_DIR}/shell"
DATA_DIR="${PROJECT_ROOT}/data"
LOG_DIR="${DATA_DIR}/logs"
SCRIPT_DIR="${DATA_DIR}/generated_scripts"

export PATH="${BIN_DIR}:${PATH}"
export CPGP_CPUS_ON_NODE="${THREADS}"
export CPGP_JOB_ID="local"

mkdir -p "${LOG_DIR}" "${SCRIPT_DIR}"

has_step() {
  local wanted="$1"
  local s
  for s in "${RUN_STEPS[@]}"; do
    [[ "${s}" == "${wanted}" ]] && return 0
  done
  return 1
}

render_template() {
  local template="$1"
  local output="$2"
  shift 2
  cp "${template}" "${output}"
  while [[ $# -gt 0 ]]; do
    local key="$1"
    local value="$2"
    sed -i "s#${key}#${value}#g" "${output}"
    shift 2
  done
  chmod +x "${output}"
}

write_combinations() {
  local project="$1"
  shift
  local -a combos=("$@")
  local pro="${DATA_DIR}/${project}"
  mkdir -p "${pro}"
  : > "${pro}/breeds_combination.txt"
  local c
  for c in "${combos[@]}"; do
    echo "${c}" >> "${pro}/breeds_combination.txt"
  done
}

run_dataset() {
  local project="$1"
  local bfile="$2"
  local phef="$3"
  local breeds="$4"
  local traits="$5"
  local trait_indices="$6"
  local all_eff="$7"
  shift 7
  local -a combos=("$@")

  local pro="${DATA_DIR}/${project}"
  mkdir -p "${pro}" "${LOG_DIR}"
  write_combinations "${project}" "${combos[@]}"

  local pi script combo combo_id
  for pi in ${trait_indices}; do
    if has_step within; then
      script="${SCRIPT_DIR}/${project}_within_${pi}.sh"
      render_template "${SHELL_DIR}/template_within.sh" "${script}" \
        "%root%" "${PROJECT_ROOT}" \
        "%project%" "${project}" \
        "%pi%" "${pi}" \
        "%ncpus%" "${THREADS}" \
        "%mem%" "25" \
        "%subdir%/" "" \
        "%phef%" "${phef}" \
        "%genf%" "${bfile}" \
        "%breeds%" "${breeds}" \
        "%traits%" "${traits}" \
        "%all_eff%" "${all_eff}" \
        "%rep%" "${REP}" \
        "%fold%" "${FOLD}" \
        "%seed%" "$(cat "${pro}/seed.txt")"
      bash "${script}"
    fi

    if has_step joint; then
      for combo in "${combos[@]}"; do
        combo_id="${combo// /_}"
        for mt in blend union; do
          script="${SCRIPT_DIR}/${project}_${mt}_${combo_id}_${pi}.sh"
          render_template "${SHELL_DIR}/template_blend_union.sh" "${script}" \
            "%root%" "${PROJECT_ROOT}" \
            "%project%" "${project}" \
            "%mt%" "${mt}" \
            "%pi%" "${pi}" \
            "%ncpus%" "${THREADS}" \
            "%mem%" "25" \
            "%subdir%/" "" \
            "%breeds%" "${combo}" \
            "%traits%" "${traits}" \
            "%GmatM%" "single" \
            "%bfileM%" ""
          bash "${script}"
        done
      done
    fi

    if has_step bayes; then
      for combo in "${combos[@]}"; do
        combo_id="${combo// /_}"
        script="${SCRIPT_DIR}/${project}_multi_${combo_id}_${pi}.sh"
        render_template "${SHELL_DIR}/template_multi.sh" "${script}" \
          "%root%" "${PROJECT_ROOT}" \
          "%project%" "${project}" \
          "%mt%" "multi" \
          "%pi%" "${pi}" \
          "%ncpus%" "${THREADS}" \
          "%mem%" "30" \
          "%subdir%/" "" \
          "%breeds%" "${combo}" \
          "%traits%" "${traits}" \
          "%bin%" "cubic" \
          "%soft%" "C" \
          "%iter%" "${MCMC_CYCLES}" \
          "%burnin%" "${BURNIN}" \
          "%seed%" "$(cat "${pro}/seed.txt")" \
          "%rg%" "0.12" \
          "%priorVar%" "pheno" \
          "%ref%" "B" \
          "%dirPre%" "" \
          "%binf%" "null" \
          "%bfileM%" ""
        bash "${script}"
      done
    fi

    if has_step mgrm; then
      for combo in "${combos[@]}"; do
        combo_id="${combo// /_}"
        for mt in blend union; do
          script="${SCRIPT_DIR}/${project}_${mt}_MGRM_${combo_id}_${pi}.sh"
          render_template "${SHELL_DIR}/template_blend_union.sh" "${script}" \
            "%root%" "${PROJECT_ROOT}" \
            "%project%" "${project}" \
            "%mt%" "${mt}" \
            "%pi%" "${pi}" \
            "%ncpus%" "${THREADS}" \
            "%mem%" "25" \
            "%subdir%/" "" \
            "%breeds%" "${combo}" \
            "%traits%" "${traits}" \
            "%GmatM%" "MTG2" \
            "%bfileM%" ""
          bash "${script}"
        done
      done
    fi
  done
}

for dataset in "${DATASETS[@]}"; do
  case "${dataset}" in
    Xie2021)
      run_dataset \
        "Xie2021" \
        "${PROJECT_ROOT}/data/Xie2021/Genotype.id.qc" \
        "${PROJECT_ROOT}/data/Xie2021/phenotypes_dmu.txt" \
        "YY LL LY" \
        "PFAI MS" \
        "1 2" \
        "'3 2 1'" \
        "YY LY" "LL LY" "YY LL LY"
      ;;
    Lee2019)
      run_dataset \
        "Lee2019" \
        "${PROJECT_ROOT}/data/Lee2019/Lee2019q" \
        "${PROJECT_ROOT}/data/Lee2019/phenotype.txt" \
        "AAN LF LIM" \
        "BWT CE CW DOC MARB MCE MILK REA SCRO STAY WWT YG YWT" \
        "1 11" \
        "'2 1'" \
        "AAN LF" "LF LIM" "AAN LF LIM"
      ;;
    *)
      echo "Unknown dataset: ${dataset}"
      exit 1
      ;;
  esac
done

echo "CPGP workflow finished. Outputs are under ${DATA_DIR}/<dataset>/."

