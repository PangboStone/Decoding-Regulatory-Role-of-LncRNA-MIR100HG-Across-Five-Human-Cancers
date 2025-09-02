import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

def Correlation(methylation_file):
    """
    Computes correlation between CpG methylation and MIR100HG expression.
    :param methylation_file: Path to the methylation file with MIR100HG expression added.
    """
    # 读取数据
    methylation_df = pd.read_csv(methylation_file, index_col=0)

    # 提取 MIR100HG 表达数据
    mir100hg_expression = methylation_df.loc["MIR100HG"]
    methylation_df = methylation_df.drop(index="MIR100HG")  # 移除 MIR100HG，剩下的就是 CpG 位点

    # 计算 Pearson 相关性
    correlations = {}
    for cpg in methylation_df.index:
        corr, p_val = stats.pearsonr(methylation_df.loc[cpg].dropna(), mir100hg_expression.dropna())
        correlations[cpg] = (corr, p_val)

    # 转换为 DataFrame
    corr_df = pd.DataFrame.from_dict(correlations, orient="index", columns=["Correlation", "P-value"])
    corr_df = corr_df.sort_values("Correlation", ascending=False)

    print(corr_df)

    # **绘制热图**
    plt.figure(figsize=(8, 5))
    sns.heatmap(corr_df[["Correlation"]], annot=True, cmap="coolwarm")
    plt.title("CpG Methylation vs. MIR100HG Expression Correlation")
    plt.show()


# 示例调用
Correlation("NTU_DATA/MIR100HG_Methylation_SKCM.csv")

