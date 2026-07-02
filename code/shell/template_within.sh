#!/bin/bash

unset R_HOME


## 交叉验证
fold=%fold%
rep=%rep%

## 路径
root=%root%
pro=${root}/data/%project%
code=${root}/code
logp=${root}/data/logs

## 表型/基因型文件
phef=%phef%
genf=%genf%

## 脚本
accur_within=${code}/shell/accur_dmu_SS_GBLUP.sh
phe_adj=${code}/shell/dmu_get_pheno_adj.sh

## 表型名
traits=(%traits%)

## 效应设置（只在品种内评估设置，多品种评估与品种内预测相比只多了品种效应）
all_eff=%all_eff%

## 随机数种子
seed=%seed%

## 品种
breeds=(%breeds%)

## 子文件夹
subdir="%subdir%/"

## 修改工作文件夹
mkdir -p ${pro}
cd ${pro} || exit

## 并行作业数
thread=${thread:-$((rep * fold))}

## 群体内评估
for pi in %pi%; do # pi=1;b=${breeds[0]}
  phedir=${pro}/${subdir}${traits[$((pi - 1))]}

  mkdir -p ${phedir}
  cd "${phedir}" || exit

  [[ -s ${phedir}/phe_adj_PBLUP.txt ]] && rm ${phedir}/phe_adj_PBLUP.txt

  for b in "${breeds[@]}"; do
    mkdir -p ${phedir}/${b}
    cd ${phedir}/${b} || exit

    ## 检查plink中有无家系id
    if [[ ! -s ${genf}.fam ]]; then
      echo "${genf}.fam not found! "
      exit 1
    elif [[ $(grep -c "${b}" ${genf}.fam) -eq 0 ]]; then
      echo "no family id ${b} in ${genf}.fam file! "
      exit 1
    fi

    ## 提取指定品种基因型(可能需要对基因型文件进行修改，所以也是复制基因型)
    echo "${b}" >tmp_fid.txt
    plink --bfile ${genf} --keep-fam tmp_fid.txt --chr-set 30 --make-bed --out ${b}
    rm tmp_fid.txt

    ## 品种相应文件
    # grep ${b} ${phef} | awk '{$NF="";print}' > ${b}_dmu.txt
    awk 'FNR==NR{a[$2];next} $1 in a' ${b}.fam ${phef} > ${b}_dmu.txt

    ## 计算校正表型
    $phe_adj \
      --phereal ${pi} \
      --bfile ${b} \
      --DIR phe_adj_PBLUP \
      --phef ${b}_dmu.txt \
      --all_eff "${all_eff}" \
      --ran_eff "1" \
      --out ${phedir}/phe_adj_PBLUP.txt \
      --append

    ## dmu计算准确性
    $accur_within \
      --label ${b} \
      --phef ${b}_dmu.txt \
      --DIR within \
      --bfile ${b} \
      --all_eff "${all_eff}" \
      --ran_eff "1" \
      --seed ${seed} \
      --rep ${rep} \
      --fold ${fold} \
      --phereal ${pi} \
      --tbvf ${phedir}/phe_adj_PBLUP.txt \
      --thread ${thread} \
      --out accur_GBLUP.txt &>${logp}/real_${traits[$((pi - 1))]}_${b}.log
  done
done

echo "job finished! "


