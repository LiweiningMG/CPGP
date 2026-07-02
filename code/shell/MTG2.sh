#!/usr/bin/env bash

## 鑷姩鍖栬绠楀鍝佺/缇や綋鍩哄洜缁勫叧绯荤煩闃碉紙GRM锛夊苟澶勭悊 NaN 闂
## 鏍稿績閫昏緫锛?
## 1. 鑷姩鏍￠獙 .fam 绗竴鍒楋紝鎻愬彇 FID 骞剁敓鎴?.pop 鍝佺鍒嗙粍鏂囦欢
## 2. 鑷姩鐢熸垚 MTG2 鎵€闇€鐨?.rtmx 鍙傛暟鍗?
## 3. 鎸夌兢浣撴媶鍒嗗苟瀵绘壘鍗曟€?SNP -> 鍓旈櫎鍚堝苟鏁版嵁涓殑鍗曟€?SNP -> 杩愯 MTG2
## 闇€瑕佺敤鍒扮殑杞欢锛歊銆丳link1.9銆丮TG2

## 鐢ㄦ硶锛?/calc_multi_grm.sh -B merge_file -O output_prefix
# -B | --bfile    杈撳叆鐨勫寘鍚簡澶氬搧绉嶇殑鍚堝苟 PLINK 浜岃繘鍒舵枃浠跺墠缂€
# -O | --out      杈撳嚭鐨?GRM 鏂囦欢鍓嶇紑
# -C | --chr-set  鏌撹壊浣撴暟鐩紝榛樿涓?29 (鐗?缇?锛岀尓鍙涓?30
# --scale         鍚勫搧绉嶇殑灏哄害鍥犲瓙锛岀┖鏍煎垎闅斿苟鍔犲弻寮曞彿锛屽 "0 -0.75" (榛樿鍏ㄤ负 0)
# --var-calc      鍩哄洜棰戠巼鏂瑰樊璁＄畻鏂瑰紡 (榛樿 2锛屽嵆 2p(1-p))
# --keep-tmp      淇濈暀涓棿鐢熸垚鐨勪复鏃舵枃浠讹紙榛樿娓呯悊锛?
# -h | --help     甯姪

###################  鍙傛暟澶勭悊  #####################
####################################################
TEMP=$(getopt -o h,B:,O:,C: --long bfile:,out:,chr-set:,scale:,var-calc:,code_path:,keep-tmp,inv,help \
              -n 'calc_multi_grm' -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"

while true; do
  case "$1" in
    -B | --bfile )   bfile="$2";     shift 2 ;; ## PLINK 鏂囦欢鍓嶇紑
    -O | --out )     out="$2" ;      shift 2 ;; ## 杈撳嚭鍓嶇紑
    -C | --chr-set ) chr_set="$2";   shift 2 ;; ## 鏌撹壊浣撹缃?[30]
    --code_path )    code_path="$2"; shift 2 ;; ## 鍏朵粬鑴氭湰璺緞
    --scale )        scale_str="$2"; shift 2 ;; ## 鏇存敼灏哄害鍥犲瓙
    --var-calc )     var_calc="$2";  shift 2 ;; ## var璁＄畻鏂瑰紡
    --keep-tmp )     keep_tmp=true;  shift   ;; ## 淇濈暀涓棿鏂囦欢
    --inv )          do_inv=true;    shift   ;; ## 杈撳嚭閫嗙煩闃?
    -h | --help )    
      echo "Usage: $0 -B <bfile_prefix> -O <output_prefix> [-C 29] [--scale \"0 0\"] [--var-calc 2] [--keep-tmp]"
      exit 1 ;;
    -- ) shift; break ;;
    * ) shift; break ;;
  esac
done

## 妫€鏌ュ繀瑕佸弬鏁?
if [[ -z ${bfile} || -z ${out} ]]; then
    echo "Error: Missing required arguments."
    echo "Please use '$0 -h' to view instructions."
    exit 1
fi

## 榛樿鍙傛暟
chr_set=${chr_set:=30}
var_calc=${var_calc:=2}
scale_str=${scale_str:=""}
keep_tmp=${keep_tmp:=false}

## 杈呭姪鑴氭湰
RSCRIPT_INV="${code_path}/R/G_inv.R"
update_grm_ids="${code_path}/R/update_grm_ids.R"

## 妫€鏌ユ枃浠舵槸鍚﹀瓨鍦?
if [[ ! -f ${bfile}.fam || ! -f ${bfile}.bim || ! -f ${bfile}.bed ]]; then
    echo "Error: PLINK files for ${bfile} not found!"
    exit 1
fi

###################   鐜妫€鏌?  #####################
####################################################

## 妫€鏌?plink 鏄惁鍙墽琛?
if ! command -v plink &> /dev/null; then
    echo "Error: plink is not executable. Please run initialize.sh or add code/bin to PATH."
    exit 1
fi

## 宸ヤ綔璺緞璁剧疆
if [[ -n "${bfile}" ]]; then
    # 鎻愬彇鐩綍閮ㄥ垎
    bdir=$(dirname "${bfile}")
    # 鎻愬彇鏂囦欢鍚嶉儴鍒?
    bbase=$(basename "${bfile}")
    
    # 濡傛灉璺緞涓嶆槸 "." (褰撳墠鐩綍)锛屽垯灏濊瘯鍒囨崲
    if [[ "${bdir}" != "." ]]; then
        if cd "${bdir}"; then
            echo "Working directory changed to: ${bdir}"
            bfile="${bbase}"
        else
            echo "Error: Cannot change directory to ${bdir}"
            exit 1
        fi
    fi
fi

## 妫€鏌?mtg2 鏄惁鍙墽琛?(鍚屾牱涓ヨ皑澶勭悊)
if ! command -v mtg2 &> /dev/null; then
    echo "Error: mtg2 is not found in PATH. Please ensure it is installed or loaded."
    exit 1
fi

echo "=========================================================="
echo "          Multi-Breed GRM Calculation Pipeline            "
echo "=========================================================="

###################  1. 鏍￠獙涓庣敓鎴?.pop 鍒嗙粍鏂囦欢  ################
echo "[Step 1]: Validating fam file and generating .pop file..."

## 妫€鏌ョ涓€鍒?FID)鏄惁鍜岀浜屽垪(IID)瀹屽叏鐩稿悓
total_count=$(wc -l < ${bfile}.fam)
same_count=$(awk '$1==$2 {c++} END {print c+0}' ${bfile}.fam)

if [[ ${same_count} -eq ${total_count} ]]; then
    echo "Error: FID (Column 1) is identical to IID (Column 2) for all individuals."
    echo "       Cannot identify distinct breeds/populations. Please fix the .fam file."
    exit 1
fi

## 鎻愬彇鏈夊簭鐨勪笉閲嶅 FID
mapfile -t fids < <(awk '{print $1}' "${bfile}.fam" | sort | uniq)
num_fids=${#fids[@]}
echo "Found ${num_fids} populations (FIDs): ${fids[*]}"

## 鐢熸垚 .pop 鏂囦欢 (FID, IID, 鏁存暟鍝佺缂栧彿)
pop_file="${out}.pop"
awk -v fids_str="${fids[*]}" '
BEGIN {
    n = split(fids_str, arr, " ");
    for(i=1; i<=n; i++) map[arr[i]] = i;
}
{
    print $1, $2, map[$1]
}' ${bfile}.fam > ${pop_file}
echo "Generated population assignment file: ${pop_file}"

###################  2. 鐢熸垚 .rtmx 鍙傛暟鍗? ######################
echo ""
echo "[Step 2]: Generating .rtmx parameter file..."

## 澶勭悊灏哄害鍥犲瓙
if [[ -z "${scale_str}" ]]; then
    # 濡傛灉鏈彁渚涳紝榛樿涓哄叏 0
    scales=()
    for ((i=0; i<num_fids; i++)); do scales+=("0"); done
else
    # 灏嗚緭鍏ョ殑瀛楃涓叉寜绌烘牸鎵撴暎鎴愭暟缁?
    read -r -a scales <<< "${scale_str}"
    if [[ ${#scales[@]} -ne ${num_fids} ]]; then
        echo "Error: The number of scale factors (${#scales[@]}) does not match the number of populations (${num_fids})."
        echo "Please provide exactly ${num_fids} values, e.g., --scale \"0 -0.75\""
        exit 1
    fi
fi

rtmx_file="${out}.rtmx"
echo "${num_fids}" > ${rtmx_file}
echo "${pop_file}" >> ${rtmx_file}
for s in "${scales[@]}"; do
    echo "${s}" >> ${rtmx_file}
done
echo "${var_calc}" >> ${rtmx_file}

echo "Generated parameter file: ${rtmx_file}"
echo "---------------------------------"
cat ${rtmx_file}
echo "---------------------------------"

###################  3. 鎸夊搧绉嶆媶鍒嗗苟妫€鏌ュ崟鎬?SNP  ################
echo ""
echo "[Step 3]: Splitting individuals by FID and finding monomorphic SNPs..."

mono_snp_file="${out}_all_monomorphic.snps"
:> ${mono_snp_file}

for fid in "${fids[@]}"; do
  echo "  -> Processing Population: ${fid}"
  
  awk -v id="${fid}" '$1==id {print $1,$2}' ${bfile}.fam > tmp_${fid}.keep
  
  plink --bfile ${bfile} --keep tmp_${fid}.keep --chr-set ${chr_set} --make-bed --out tmp_${bfile}_${fid} > /dev/null 2>&1
  plink --bfile tmp_${bfile}_${fid} --chr-set ${chr_set} --freq --out tmp_${bfile}_${fid} > /dev/null 2>&1
  
  mono_count=$(awk 'NR>1 && ($5+0)==0 {c++} END {print c+0}' tmp_${bfile}_${fid}.frq)
  echo "     - MAF=0 SNPs in ${fid}: ${mono_count}"
  
  awk 'NR>1 && ($5+0)==0 {print $2}' tmp_${bfile}_${fid}.frq >> ${mono_snp_file}
done

###################  4. 姹囨€诲苟鍓旈櫎鎵€鏈夊崟鎬?SNP  ###################
echo ""
echo "[Step 4]: Consolidating and removing monomorphic SNPs..."

uniq_mono_file="${out}_uniq_monomorphic.snps"
sort ${mono_snp_file} | uniq > ${uniq_mono_file}
total_mono=$(wc -l < ${uniq_mono_file})

echo "Total unique monomorphic SNPs across ALL populations to remove: ${total_mono}"

clean_bfile="${out}_clean"
if [[ ${total_mono} -gt 0 ]]; then
  plink --bfile ${bfile} --exclude ${uniq_mono_file} --chr-set ${chr_set} --make-bed --out ${clean_bfile} > plink_clean.log 2>&1
  echo "Cleaned PLINK dataset generated: ${clean_bfile}"
else
  echo "No monomorphic SNPs found. Using original bfile."
  clean_bfile=${bfile}
fi

###################  5. 杩愯 MTG2 骞舵鏌?NaN  ####################
echo ""
echo "[Step 5]: Running MTG2 and checking for NaN values..."

grm_out="${out}.grm"
rm -f ${grm_out}

mtg2 -plink ${clean_bfile} -rtmx2 ${rtmx_file} -out ${grm_out} > mtg2.log 2>&1

if [[ ! -f ${grm_out} ]]; then
  echo "Error: MTG2 failed to generate ${grm_out}. Check mtg2.log for details."
  exit 1
fi

nan_count=$(awk '$3=="NaN" || $3=="nan" || $3=="NAN" {c++} END {print c+0}' ${grm_out})
echo "NaN rows in final GRM: ${nan_count}"

if [[ ${nan_count} -gt 0 ]]; then
  echo "Warning: There are still NaN values in the GRM!"
  echo "First 5 NaN rows:"
  awk '$3=="NaN" || $3=="nan" || $3=="NAN" {print $0}' ${grm_out} | head -n 5
else
  echo "Success: No NaN values found in the GRM."
fi

###################   Step 4.9锛氭浛鎹?GRM 涓殑 ID   ####################
echo ""
echo "[Step 4.9]: Replacing GRM indices with FAM IDs using getopt script..."
# 鍋囪 ${fam_path} 鏄綘鐨?fam 鏂囦欢璺緞
# ${grm_out} 鏄師濮嬬敓鎴愮殑 grm 鏂囦欢
updated_grm="${out}.named.grm"
# 璋冪敤涓婇潰鐨?R 鑴氭湰
$update_grm_ids --fam "${bfile}.fam" --grm "${grm_out}" --out "${grm_out}"

if [[ $? -eq 0 ]]; then
    echo "Success: IDs replaced."
else
    echo "Error: ID replacement failed."
    exit 1
fi

###################   Step 5 寤朵几锛欸RM 姹傞€?  ####################
if [[ ${do_inv} == true ]]; then
    echo ""
    echo "[Step 5.1]: Inverting the GRM using R..."

    if [[ ! -f ${RSCRIPT_INV} ]]; then
        echo "Error: Inversion R script not found at ${RSCRIPT_INV}"
        exit 1
    fi

    # 璋冪敤 R 鑴氭湰
    # 浼犲叆鍙傛暟锛氳緭鍏ユ枃浠讹紝杈撳嚭鍓嶇紑
    "${RSCRIPT_INV}" --grm "${grm_out}" --out "${out}" --tol 1e-6
    
    if [[ $? -eq 0 ]]; then
        echo "Success: Inverse GRM generated as ${out}.grm.inv"
    else
        echo "Error: R inversion failed."
        exit 1
    fi
fi

###################  6. 娓呯悊涓棿鏂囦欢  ###########################
echo ""
echo "[Step 6]: Cleanup..."

if [[ ${keep_tmp} == false ]]; then
  rm -f tmp_*
  rm -f ${mono_snp_file} ${uniq_mono_file}
  if [[ ${total_mono} -gt 0 ]]; then
    rm -f ${clean_bfile}.bed ${clean_bfile}.bim ${clean_bfile}.fam plink_clean.log
  fi
  echo "Temporary PLINK and SNP files removed."
  echo "Kept configuration files: ${pop_file}, ${rtmx_file}"
else
  echo "Temporary files kept as requested."
fi

echo "==================== Pipeline Finished ===================="


