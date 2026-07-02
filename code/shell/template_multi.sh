#!/usr/bin/env bash

unset R_HOME

root=%root%
pro=${root}/data/%project%
logp=${root}/data/logs
code=${root}/code

bfileM="%bfileM%"
bfileM_arg=()
[[ -n "${bfileM}" ]] && bfileM_arg=(--bfileM "${bfileM}")

breeds=(%breeds%)
traits=(%traits%)
nsnp_win=100
seed=%seed%
subdir="%subdir%/"
accur_multi=${code}/shell/accur_multi_breed.sh

mkdir -p "${logp}"
cd "${pro}" || exit 1

for pi in %pi%; do
  ti=${traits[$((pi - 1))]}
  phep="${pro}/${subdir}${ti}"
  cd "${phep}" || exit 1

  suffix=$(IFS=_ ; echo "_${breeds[*]}")
  if [[ -d ${phep}/multi${suffix} ]]; then
    echo "${phep}/multi${suffix} exists; remove it before rerunning."
    exit 0
  fi

  if [[ -f random_seed.txt ]]; then
    seed=$(cat random_seed.txt)
  elif [[ -z "${seed}" ]]; then
    seed=$RANDOM
    echo "${seed}" > random_seed.txt
  fi

  rep=$(find ./${breeds[0]}/val1/rep* -type d 2>/dev/null | wc -l)
  fold=$(find ./${breeds[0]}/val*/rep1 -type d 2>/dev/null | wc -l)
  thread=${thread:-$((rep * fold))}

  for soft in %soft%; do
    for bin in %bin%; do
      "${accur_multi}" \
        --pops "${breeds[*]}" \
        "${bfileM_arg[@]}" \
        --tbvf "${phep}/phe_adj_PBLUP.txt" \
        --fold "${fold}" \
        --rep "${rep}" \
        --type multi \
        --software "${soft}" \
        --bin "${bin}" \
        --phereal "${pi}" \
        --nsnp_win "${nsnp_win}" \
        --thread "${thread}" \
        --rg %rg% \
        --priorVar %priorVar% \
        --binf %binf% \
        --iter %iter% \
        --burnin %burnin% \
        --ref %ref% \
        --dirPre %dirPre% \
        --res_const \
        --rg_local \
        --suffix \
        --seed "${seed}" \
        --code "${code}" \
        --out accur_bayes >"${logp}/%project%_${ti}_multi${suffix}_${soft}_${bin}_${CPGP_JOB_ID:-local}.log" 2>&1
    done
  done
done

echo "finish normally"
