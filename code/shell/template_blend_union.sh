#!/usr/bin/env bash

unset R_HOME

root=%root%
pro=${root}/data/%project%
logp=${root}/data/logs
code=${root}/code

breeds=(%breeds%)
traits=(%traits%)

bfileM="%bfileM%"
bfileM_arg=()
[[ -n "${bfileM}" ]] && bfileM_arg=(--bfileM "${bfileM}")

subdir="%subdir%/"
accur_multi=${code}/shell/accur_multi_breed.sh

mkdir -p "${logp}"
cd "${pro}" || exit 1

for pi in %pi%; do
  ti=${traits[$((pi - 1))]}
  phep="${pro}/${subdir}${ti}"
  cd "${phep}" || exit 1

  rep=$(find ./${breeds[0]}/val1/rep* -type d 2>/dev/null | wc -l)
  fold=$(find ./${breeds[0]}/val*/rep1 -type d 2>/dev/null | wc -l)
  thread=${thread:-$((rep * fold))}

  for m in %mt%; do
    suffix=$(IFS=_ ; echo "_${breeds[*]}")
    scenario="${m}${suffix}"
    if [[ "%GmatM%" == "single" ]]; then
      scenario="${m}_SGRM${suffix}"
    elif [[ "%GmatM%" == "MTG2" ]]; then
      scenario="${m}_MGRM${suffix}"
    fi
    if [[ -d ${phep}/${scenario} ]]; then
      echo "${phep}/${scenario} exists; remove it before rerunning."
      exit 0
    fi

    "${accur_multi}" \
      --pops "${breeds[*]}" \
      "${bfileM_arg[@]}" \
      --type "${m}" \
      --phereal "${pi}" \
      --GmatM %GmatM% \
      --fold "${fold}" \
      --rep "${rep}" \
      --suffix \
      --thread "${thread}" \
      --tbvf "${phep}/phe_adj_PBLUP.txt" \
      --code "${code}" \
      --out accur_GBLUP >"${logp}/%project%_${ti}_${scenario}_${CPGP_JOB_ID:-local}.log" 2>&1
  done
done

echo "finish normally"
