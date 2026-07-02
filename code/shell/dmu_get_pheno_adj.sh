#!/usr/bin/bash

## date
## liwn: 2022-07-25
## 用于获取校正所有固定效应和非遗传效应的校正表型，用于准确性计算(Christensen et al., 2012)

################  命令行参数处理  ##################
####################################################
## NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
## 参数名
TEMP=$(getopt -o 4vh\?p:b:m: --long append,ped_var,dmu4,help,phef:,bfile:,gmat:,gidf:,pedf:,all_eff:,ran_eff:,add_rf:,invA:,varf:,phereal:,miss:,num_int:,code:,DIR:,alpha:,out:,intercept,debug \
    -n 'javawrap' -- "$@")
if [ $? != 0 ]; then
    echo "Open script $0 to view instructions"
    echo "Terminating..." >&2
    exit 1
fi
eval set -- "$TEMP"

## 解析参数
while true; do
    case "$1" in
    -p | --phef )     phef="$2";   shift 2 ;;  # 表型文件
    -b | --bfile )    bfile="$2";  shift 2 ;;  # plink二进制文件前缀
        --pedf )     pedf="$2";    shift 2 ;;  # 系谱文件(不提供则运行GBLUP)
    -m | --gmat )     gmat="$2";   shift 2 ;;  # 用户提供关系矩阵或逆矩阵(id id value)
        --gidf )     gidf="$2";    shift 2 ;;  # 基因型个体id，与用户指定的G阵文件中个体id一致
        --all_eff )  all_eff="$2"; shift 2 ;;  # DIR中$MODEL第3行(所有效应)，前3位不用，只需所有效应所在的列，如"2 3 1"
        --ran_eff )  ran_eff="$2"; shift 2 ;;  # DIR中$MODEL第4行(随机效应分)，第1位不用，只需所有随机效应所在分组，如"1"
        --add_rf )   add_rf="$2";  shift 2 ;;  # 加性效应所在分组
        --add_sol )  add_sol="$2"; shift 2 ;;  # 加性效应在SOL文件第一列中的编号
        --invA )     invA="$2";    shift 2 ;;  # A逆构建方式(1/2/3/4/6)，1为考虑近交，2为不考虑近交，其他见DMU说明书
        --varf )     varf="$2";    shift 2 ;;  # 方差组分文件(如没提供则用系谱/SNP信息估计)
        --phereal )  phereal="$2"; shift 2 ;;  # 表型在表型文件中实数列的位置
        --miss )     miss="$2";    shift 2 ;;  # 缺失表型表示符
        --num_int )  num_int="$2"; shift 2 ;;  # 整型列列数
        --code )     code="$2";    shift 2 ;;  # 代码路径
        --DIR )      DIR="$2";     shift 2 ;;  # 参数卡文件前缀
        --alpha )    alpha="$2";   shift 2 ;;  # G阵校正系数(是否考虑近交 )
        --out )      out="$2";     shift 2 ;;  # 输出校正表型文件名
        --append )   append=true;  shift   ;;  # 在结果文件中追加，而不是覆盖
        --debug)     debug=true;   shift   ;;  # 在整数列最后一列后增加一列群体均值
        --intercept) mean=true;    shift   ;;  # 在整数列最后一列后增加一列群体均值
    -v | --ped_var )  ped_var=true; shift   ;;  # 用系谱信息估计方差组分
    -4 | --dmu4 )     dmu4=true;    shift   ;;  # 用DMU4模型估计育种值
    -h | --help | -\? )  echo "Open script $0 to view instructions" && exit 1 ;;
    -- ) shift; break ;;
    * ) break ;;
    esac
done

## 参数默认值
ran_eff=${ran_eff:=1}                      ## 默认只有加性遗传一个效应
all_eff=${all_eff:="2 1"}                  ## 默认只有群体均值一个固定效应，且在表型文件第二列(全为1)
phereal=${phereal:=1}                      ## 表型列
add_rf=${add_rf:=1}                        ## 加性随机效应所在组
miss=${miss:=-99}                          ## 缺失表型表示
invA=${invA:=1}                            ## A逆构建方式(是否考虑近交)
alpha=${alpha:=0.05}                       ## G阵校正系数(是否考虑近交)
DIR=${DIR:=phe_adj}                        ## 参数卡文件前缀
out=${out:=phe_adj.txt}                    ## 输出文件名
append=${append:=false}                    ## 结果附加在已有文件上，而不是覆盖
varf=${varf:=}                             ## 避免vscode报错

####################
## 需要调用的软件
# R plink gmatrix

## 避免执行R脚本时的警告("ignoring environment value of R_HOME")
unset R_HOME

## 检查必要参数是否提供
if [[ ! -s ${phef} ]]; then
    echo "phenotype file ${phef} not found! "
    exit 1
elif [[ ! -s ${bfile}.fam && ! -s ${gidf} && ! -s ${pedf} ]]; then
    echo "plink file ${bfile}.fam, pedigree file ${pedf} or genotyped individuals id file ${gidf} not found! "
    exit 1
elif [[ -s ${gmat} && ! -s ${gidf} ]]; then
    echo "genotyped individuals id file ${gidf} not found! "
    exit 1
fi

## 主文件夹
workdir=$(pwd)

## 脚本
keep_phe_gid=${code}/R/pheno_all_genoid.R
miss_phe=${code}/R/pheno_miss.R
phe_adj=${code}/R/pheno_adj.R

## 检查文件中是否含有字符
check_alphabet() {
    [[ ! -s ${1} ]] && echo "${1} not found! " && exit 1
    if [[ ${2} ]]; then
        ## 科学计数法表示的数跳过("e-")
        NotNum=$(awk -vl=${2} '{print $l}' ${1} | grep -v "[0-9]e\-[0-9]" | grep -c "[a-zA-Z]")
    else
        NotNum=$(cat ${1} | grep -v "[0-9]e\-[0-9]" | grep -c "[a-zA-Z]")
    fi

    if [[ ${NotNum} -gt 0 ]]; then
        echo "Non numeric characters exist in ${1} file, please check! "
        exit 1
    fi
}

## 表型文件整型、实型变量列数
ncol=$(awk 'END{print NF}' ${phef})
for i in $(seq 1 ${ncol}); do
    dot=$(awk -vl=${i} '{print $l}' ${phef} | grep -c "\.")
    [[ ${dot} -gt 0 ]] && num_int=$((i - 1)) && break
done
num_real=$(($(awk 'END{print NF}' ${phef}) - num_int))

## 效应个数
nA=$(echo ${all_eff} | awk '{print NF}')
nR=$(echo ${ran_eff} | awk '{print NF}')

###################  表型文件处理  #####################
########################################################
## 剔除缺失表型个体
echo "remove individuals missing phenotypes in the phenotype file"
[[ -s ${bfile}.fam ]] && option="--map ${bfile}.fam"
$miss_phe \
    --file ${phef} \
    --col $((phereal + num_int)) \
    --miss ${miss} \
    --missid ${workdir}/miss_phe.id \
    --out "${workdir}/pheno_adj.txt" \
    ${option}

phef=${workdir}/pheno_adj.txt
if [[ ${bfile} ]]; then
    ## 在基因型文件中剔除缺失表型个体
    if [[ -s ${workdir}/miss_phe.id ]]; then
        n_miss_phe=$(cat ${workdir}/miss_phe.id | wc -l)
        echo "remove ${n_miss_phe} individuals with the missing value in plink files"
        plink --bfile ${bfile} --chr-set 30 --make-bed --out ${bfile}.org >${workdir}/plink.log
        plink --bfile ${bfile}.org --chr-set 30 --remove ${workdir}/miss_phe.id --make-bed --out ${bfile} >>${workdir}/plink.log
    fi
fi
## 截距项列设置(整列设为"1")
if [[ ${mean} ]]; then
    ## 效应参数
    ((nA++))
    ((num_int++))
    all_eff="${num_int} ${all_eff}"

    awk -v column="${num_int}" -v value="1" '
    BEGIN {
        FS = OFS = " ";
    }
    {
        for ( i = NF + 1; i > column; i-- ) {
            $i = $(i-1);
        }
        $i = value;
        print $0;
    }  
    ' ${phef} >${phef}.tmp
    echo "add populations mean column ${num_int} in the phenotype file."
    mv ${phef}.tmp ${phef}
fi

###################  dmu参数卡模板  ####################
########################################################
[[ -s ${DIR}.DIR ]] && echo "warn: ${DIR}.DIR will be overwrited! "
{
    echo "\$COMMENT"
    echo "creating phenotypes corrected for fixed effects and non-genetic random effects"
    echo "\$ANALYSE %ANALYSE%"
    echo "\$DATA  ASCII (${num_int}, ${num_real}, ${miss}) %phef%"
    echo -e "\$MODEL\n1\n0\n${phereal} 0 ${nA} ${all_eff}\n${nR} ${ran_eff}\n0\n0"
    echo -e "\$VAR_STR %VAR_STR%"
    echo -e "\$PRIOR %PRIOR%"
    echo -e "\$RESIDUALS ASCII"
    echo -e "\$SOLUTION"
} >${DIR}.DIR

##################  文件合法性检查  ####################
########################################################
## 表型文件(是否有非数字字符)
sed -i "s/na/${miss}/Ig" ${phef}
check_alphabet ${phef}
[[ -s ${pedf} ]] && check_alphabet ${pedf}
[[ -s ${gmat} ]] && check_alphabet ${gmat}
[[ -s ${bfile}.fam ]] && check_alphabet ${bfile}.fam 2

###################  方差组分估计  #####################
########################################################
## 检查是否需要估计方差组分(用系谱)
if [[ ${ped_var} ]]; then
    ## 检查文件
    [[ ! -s ${pedf} ]] && echo "${pedf} not found! " && exit 1

    ## 替换参数卡信息
    sed 's#%ANALYSE%#1 1 0 0#g' ${DIR}.DIR >ped_var.DIR
    sed -i '/\$RESIDUALS.*/d' ped_var.DIR
    sed -i '/\$SOLUTION.*/d' ped_var.DIR
    sed -i "s#%VAR_STR%#${add_rf} PED ${invA} ASCII ${pedf}#g" ped_var.DIR
    ## 方差组分
    if [[ -s ${varf} ]]; then
        sed -i "s#%PRIOR%#${varf}#g" ped_var.DIR
    else
        sed -i '/\$PRIOR.*/d' ped_var.DIR
    fi

    ## 使用AIREML估计方差组分
    run_dmuai ped_var

    ## 方差组分文件
    sed -i 's#%PRIOR%#ped_var.PAROUT#g' ${DIR}.DIR
elif [[ -s ${varf} ]]; then
    sed -i 's#%PRIOR%#\${varf}#g' ${DIR}.DIR
else
    sed -i '/\$PRIOR.*/d' ${DIR}.DIR
fi

###################  方差组分结构  ####################
#######################################################
if [[ -s ${bfile}.fam || -s ${gmat} ]] && [[ ! -s ${pedf} ]]; then
    ## GBLUP
    method=GBLUP
    gmat=${gmat:=full.agiv.id_fmt} ## 3列格式G逆阵文件(id id value)
    sed -i "s#%VAR_STR%#${add_rf} GREL ASCII ${gmat}#g" ${DIR}.DIR

    ## 只保留有基因型个体的表型
    $keep_phe_gid \
        --famf "${bfile}.fam" \
        --phef ${phef} \
        --out ${phef}_gid
    phef=${phef}_gid
    add_sol=3
elif [[ ! -s ${bfile}.fam && -s ${pedf} ]]; then
    ## PBLUP
    method=PBLUP
    sed -i "s#%VAR_STR%#${add_rf} PED ${invA} ASCII ${pedf}#g" ${DIR}.DIR
    add_sol=4
else
    method=ssGBLUP
    gmat=${gmat:=full.agrm.id_fmt} ## 3列格式G阵文件(id id value)
    gidf=${gidf:=full.id}          ## 1列 id
    sed -i "s#%VAR_STR%#${add_rf} PGMIX ${invA} ASCII ${pedf} ${gidf} ${gmat} ${alpha} G-ADJUST#g" ${DIR}.DIR
    add_sol=4
fi

###################  关系矩阵构建  #####################
########################################################
## 构建G逆矩阵(同时输出基因型个体id) gmatrix软件
if [[ ! -s ${gmat} && ${bfile} ]]; then
    [[ ${method} == "GBLUP" ]] && inv=" --inv" || inv=""
    echo "Read the plink bed file and Calculate the additive G matrix..."
    gmatrix --bfile ${bfile} --grm agrm --out full ${inv} >gmatrix.log

    if [[ $? -ne 0 ]]; then
        echo "G matrix calculate error! "
        exit 1
    else
        echo "G matrix created."
    fi
fi

#################  育种值估计  ###################
##################################################
echo "estimating breeding values using the ${method} model..."
sed -i "s#%phef%#${phef}#g" ${DIR}.DIR
if [[ ! ${debug} ]]; then
    if [[ ${dmu4}  ]]; then
        [[ ! -s ${varf} ]] && echo "${varf} not found! " && exit 1
        ## 运行dmu4
        sed -i 's#%ANALYSE%#11 9 0 0#g' ${DIR}.DIR
        run_dmu4 ${DIR}
        [[ $? -ne 0 ]] && echo "error in dmu4! " && exit 1
    else
        sed -i 's#%ANALYSE%#1 1 0 0#g' ${DIR}.DIR
        [[ ! ${debug} ]] && run_dmuai ${DIR}
        [[ $? -ne 0 ]] && echo "error in dmuai! " && exit 1
    fi
fi

#################  计算校正表型  #################
##################################################
echo "calculating adjusted phenotype."
$phe_adj \
    --DIR ${DIR} \
    --phe ${phef} \
    --add_sol ${add_sol} \
    --out ${out} \
    --append ${append}
