#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"
CODE_DIR="${PROJECT_ROOT}/code"
BIN_DIR="${CODE_DIR}/bin"
DATA_DIR="${PROJECT_ROOT}/data"
LOG_DIR="${DATA_DIR}/logs"
SCRIPT_DIR="${DATA_DIR}/generated_scripts"

mkdir -p "${LOG_DIR}" "${SCRIPT_DIR}"

missing=0
check_file() {
  local f="$1"
  if [[ ! -s "${f}" ]]; then
    echo "Missing: ${f}"
    missing=1
  else
    echo "OK: ${f}"
  fi
}

check_exe() {
  local f="$1"
  if [[ ! -s "${f}" ]]; then
    echo "Missing executable: ${f}"
    missing=1
  else
    chmod +x "${f}" || true
    echo "OK: ${f}"
  fi
}

echo "Checking basic data files..."
check_file "${DATA_DIR}/Xie2021/Genotype.id.qc.bed"
check_file "${DATA_DIR}/Xie2021/Genotype.id.qc.bim"
check_file "${DATA_DIR}/Xie2021/Genotype.id.qc.fam"
check_file "${DATA_DIR}/Xie2021/phenotypes_dmu.txt"
check_file "${DATA_DIR}/Xie2021/seed.txt"
check_file "${DATA_DIR}/Lee2019/Lee2019q.bed"
check_file "${DATA_DIR}/Lee2019/Lee2019q.bim"
check_file "${DATA_DIR}/Lee2019/Lee2019q.fam"
check_file "${DATA_DIR}/Lee2019/phenotype.txt"
check_file "${DATA_DIR}/Lee2019/seed.txt"

echo "Checking executable programs..."
for exe in plink gmatrix dmu1 dmuai run_dmu4 run_dmuai mtg2 mbBayesABLD ldblock LD_mean_r2; do
  check_exe "${BIN_DIR}/${exe}"
done

echo "Checking scripts..."
for f in "${CODE_DIR}"/*.sh "${CODE_DIR}/shell"/*.sh; do
  [[ -e "${f}" ]] || continue
  chmod +x "${f}" || true
  bash -n "${f}"
  echo "OK: ${f}"
done

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Missing command: Rscript"
  missing=1
else
  echo "OK: Rscript"
fi

if [[ ${missing} -ne 0 ]]; then
  echo "Initialization finished with missing files or commands. Please fix them before running main.sh."
  exit 1
fi

echo "Initialization finished successfully. Edit PROJECT_ROOT in main.sh, then run: bash main.sh"

