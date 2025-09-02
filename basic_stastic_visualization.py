import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_mir100hg_expression_distribution(expression_file):
    """
    Plots the distribution of MIR100HG expression in LUAD samples.
    :param expression_file: Path to the TOIL_Filtered_Expression_XXXX.csv file.
    """
    # 读取基因表达数据，确保索引是字符串
    expression_df = pd.read_csv(expression_file, index_col=0)

    mir100hg_expression = expression_df.loc["MIR100HG"]

    plt.figure(figsize=(12, 5))

    # 直方图（查看分布情况）
    plt.subplot(1, 2, 1)
    sns.histplot(mir100hg_expression, bins=30, kde=True, color='skyblue')
    plt.title("MIR100HG Expression Distribution")
    plt.xlabel("Expression Level (TPM)")
    plt.ylabel("Frequency")

    # 箱线图（查看异常值）
    plt.subplot(1, 2, 2)
    sns.boxplot(x=mir100hg_expression, color='lightcoral')
    plt.title("MIR100HG Expression Boxplot")

    plt.tight_layout()
    plt.show()

    # 输出统计信息
    print(mir100hg_expression.describe())

# 示例：可视化 LUAD 数据的 MIR100HG 表达分布
expression_file = "./NTU_DATA/MIR100HG_Methylation_PRAD.csv"
plot_mir100hg_expression_distribution(expression_file)
