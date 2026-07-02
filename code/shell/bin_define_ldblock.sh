#!/usr/bin/bash
## 使用ldblock()将基因组根据ld信息划分成近似独立的ld块
## ldblock: https://github.com/cadeleeuw/lava-partitioning

###################  参数处理  #####################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o h --long bfile:,type:,spar:,diff:,maf:,win:,jobs:,minSize:,code:,out:,keepTmp,plot \
  -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
  echo "Terminating..." >&2
  exit 1
fi
eval set -- "$TEMP"

## 解析参数
while true; do
    case "$1" in
    --bfile )   bfile="$2";    shift 2 ;; ## PLINK_PREFIX
    --type )    type="$2";     shift 2 ;; ## Method for determining bins [lava/cubic]
    --spar )    spar="$2";     shift 2 ;; ## Smoothness (0-1) when fitting curves [0.2]
    --diff )    diff="$2";     shift 2 ;; ##  [0.05]
    --maf )     maf="$2";      shift 2 ;; ## MAF filtering threshold [0.01]
    --geno )    geno="$2";     shift 2 ;; ## remove SNPs with missing call rates [0.2]
    --mind )    mind="$2";     shift 2 ;; ## remove individuals with missing call rates [0.2]
    --win )     win="$2";      shift 2 ;; ## Correlation window in number of SNPs [50]
    --jobs )    jobs="$2";     shift 2 ;; ## Number of parallel jobs [5]
    --minSize ) minSize="$2";  shift 2 ;; ## The minimum size in number of SNPs a block must have [50]
    --out )     out="$2";      shift 2 ;; ## The prefix used for output files [ldblock]
    --nchr )    nchr="$2";     shift 2 ;; ## number of chromosomes [30]
    --keepTmp ) keepTmp=true;  shift   ;; ## Keep intermediate files [false]
    --plot )    plot="--plot"; shift   ;; ## plot mean r2 [FALSE]
    -h | --help ) grep " ;; ## " $0 && exit 1 ;;
    -- ) shift; break ;;
    * ) break ;;
    esac
done

## 软件加载
if [ ! "$(command -v ldblock)" ]; then
  echo "command \"ldblock\" does not exists on system! "
  exit 1
fi
if [ ! "$(command -v plink)" ]; then
  echo "command \"plink\" does not exists on system! "
  exit 1
fi
if [[ ${type} == "cubic" && ! "$(command -v LD_mean_r2)" ]]; then
  echo "command \"LD_mean_r2\" does not exists on system! "
  exit 1
fi

unset R_HOME

## 软件/脚本
code=${code:=$(cd "$(dirname "$0")/.." && pwd)}
Rcubic=${code}/R/cubic_smoothing_block.R

## 其他默认参数
out=${out:=ldblock.txt}
type=${type:=cubic}
spar=${spar:=0.2}
diff=${diff:=0.05}
maf=${maf:=0}
frq=${frq:=0}
geno=${geno:=0.2}
mind=${mind:=0.2}
win=${win:=50}
minSize=${minSize:=50}
jobs=${jobs:=5}
nchr=${nchr:=30}

## 检查文件格式
if [[ ! -s ${bfile}.bim ]]; then
  if [[ -s ${bfile}.map ]]; then
    plink --file ${bfile} --chr-set ${nchr} --make-bed --out ${bfile}
  else
    echo "plink file ${bfile}.map(.bim) not found! "
    exit 1
  fi
fi

## 并行作业数
job_pool=${code}/shell/job_pool.sh
source ${job_pool}
job_pool_init ${jobs} 0

## 随机种子，防止不同进程之间干扰
seed=$RANDOM

## 每条染色体进行区间划分
chrs=$(awk '{print $1}' ${bfile}.bim | sort -n | uniq)
nchr=$(echo ${chrs} | tr " " "\n" | wc -l)
for chr in ${chrs}; do
  ## 提取染色体信息
  plink --bfile ${bfile} --chr-set ${nchr} --chr ${chr} --make-bed --out ld_block_tmp_${chr}.${seed} >/dev/null

  ## SNP数目少于区间最少标记数，则跳过，未检测到结果文件会将这条染色体单独划分在一个区间
  [[ $(wc -l <ld_block_tmp_${chr}.${seed}.bim) -lt ${minSize} ]] && continue

  ## 区间划分
  if [[ ${type} == "lava" ]]; then
    ldblock \
      ld_block_tmp_${chr}.${seed} \
      -frq ${frq} \
      -win ${win} \
      -min-size ${minSize} \
      -out ld_block_tmp_${chr}.${seed}
  else
    LD_mean_r2 \
      --bfile ld_block_tmp_${chr}.${seed} \
      --maf ${maf} \
      --geno ${geno} \
      --mind ${mind} \
      --win ${win} \
      --min ${minSize} \
      --out ld_mean.${seed}
    [[ -f ld_mean.${seed}_chr${chr}.txt ]] && mv ld_mean.${seed}_chr${chr}.txt ld_block_tmp_${chr}.${seed}.breaks
  fi
done

## 等待后台程序运行完毕
job_pool_wait
job_pool_shutdown

## 准备区间文件
: >ld_block_tmp.${seed}.breaks
for chr in ${chrs}; do
  block_file=ld_block_tmp_${chr}.${seed}.breaks

  if [[ -s ${block_file} ]]; then
    if [[ ${type} == "lava" ]]; then
      nblock=$(($(wc -l <${block_file}) - 2))

      ## 准备区间文件
      mapfile -t -O 1 start < <(awk -v lines="$(wc -l <${block_file})" 'NR>1 && NR<lines {print $6}' ${block_file})
      mapfile -t -O 1 stop < <(awk 'NR>2 {print $6-1}' ${block_file})
      mapfile -t -O 1 nsnp < <(awk 'NR > 2 {print $5 - prev} {prev = $5}' ${block_file})

      ## 输出到文件
      for line in $(seq ${nblock}); do
        echo "${chr} ${start[line]} ${stop[line]} ${nsnp[line]}" >>ld_block_tmp.${seed}.breaks
      done
    else
      ${Rcubic} \
        --r2 ${block_file} \
        --bim ld_block_tmp_${chr}.${seed}.bim \
        --spar ${spar} \
        --diff ${diff} \
        --out ld_block_tmp.${seed}.cubic \
        ${plot}
      cat ld_block_tmp.${seed}.cubic >>ld_block_tmp.${seed}.breaks
    fi
  else
    starti=$(head -n 1 ld_block_tmp_${chr}.${seed}.bim | awk '{print $4}')
    stopi=$(tail -n 1 ld_block_tmp_${chr}.${seed}.bim | awk '{print $4}')
    echo "${chr} ${starti} ${stopi} $(wc -l <ld_block_tmp_${chr}.${seed}.bim)" >>ld_block_tmp.${seed}.breaks
  fi
done

## 增加序列行
awk '{print NR, $0}' ld_block_tmp.${seed}.breaks >${out}

## 检查标记数目
nsnp_bin=$(awk '{sum+=$5}END{print sum}' ${out})
if [[ ${nsnp_bin} -ne $(wc -l <${bfile}.bim) ]]; then
  echo "number of snp in bin file (${nsnp_bin}) are not equal to bfile ($(wc -l <${bfile}.bim))! "
  # rm ${out}
  rm ld_block_tmp*
  exit 2
fi

## 插入标题行
sed -i '1i LOC CHR START STOP nSNP' ${out}

## 删除中间文件
[[ ! ${keepTmp} ]] && rm ld_block_tmp*

## 报告
echo "number of blocks: $(sed '1d' ${out} | wc -l)"
echo "blocks information file output to: ${out}"

win=50
type=cubic
maf=-0.01
