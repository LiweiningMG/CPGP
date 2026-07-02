# 可执行程序说明

本目录用于保存复现本项目所需的可执行程序。当前流程可能调用以下程序：

```text
plink          PLINK，用于基因型处理、PCA 和 LD 计算
mtg2           MTG2，用于构建多群体基因组关系矩阵
dmuai          DMUAI，用于方差组分估计
dmu1           DMU1，用于 GBLUP 育种值预测
run_dmuai      DMUAI 辅助运行脚本
run_dmu4       DMU 辅助运行脚本
gmatrix        基因组关系矩阵相关辅助程序
ldblock        LD 区块划分辅助程序
LD_mean_r2     LD 衰减统计辅助程序
mbBayesABLD    本研究使用的多群体贝叶斯基因组预测程序
QMSim_selected 原始流程中保留的兼容性辅助程序
```

这些程序仅作为复现本项目分析结果的辅助文件提供。任何其他用途都应遵守对应软件作者、开发团队或发布机构的许可协议。若用户已经在系统环境中安装了这些程序，也可以在 `main.sh` 中修改路径，优先调用系统版本。

在运行分析前，请确保本目录下程序具有执行权限：

```bash
chmod +x code/bin/*
```

