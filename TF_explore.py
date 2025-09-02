import pandas as pd
import scipy.stats as stats
import networkx as nx
import matplotlib.pyplot as plt

tf_target_df = pd.read_csv("NTU_DATA/gene_attribute_edge_MIR100HG.csv")
tf_target_df = tf_target_df["target"]
# print(tf_target_df)



"""
✅ gene_attribute_edges.txt.gz数据集列名： ['GeneSym', 'NA', 'GeneID', 'GeneSym.1', 'NA.1', 'GeneID.1', 'weight']
文件中：
✅  0 个可能调控 MIR100HG 的 TFs
✅  57个被 MIR100HG 调控的 TFs

背景补充：
MIR100HG 是一种长链非编码 RNA（lncRNA），其主要功能是通过调控其他基因的表达来发挥作用。
研究现状：目前的研究主要集中在 MIR100HG 如何通过调控其他基因来影响细胞功能和疾病进程，而对其自身被调控的机制研究相对较少
"""


def filter_targets(expression_file, target_genes):
    """
    Filters target genes that have significant correlation with MIR100HG.
    :param expression_file: Path to gene expression data.
    :param target_genes: List of target genes regulated by MIR100HG.
    :return: List of significantly correlated target genes.
    """
    # 读取基因表达数据
    expression_df = pd.read_csv(expression_file, index_col="HGNC_symbol")
    mir100hg_expression = expression_df.loc["MIR100HG"]
    significant_targets = []

    # 计算相关性
    for gene in target_genes:
        if gene in expression_df.index:
            corr, p_val = stats.pearsonr(expression_df.loc[gene].dropna(), mir100hg_expression.dropna())
            if p_val < 0.05 and abs(corr) > 0.5: # 这里条件的判定很关键
                significant_targets.append((gene, corr, p_val))

    # 转换为 DataFrame
    significant_df = pd.DataFrame(significant_targets, columns=["Gene", "Correlation", "P-value"])
    significant_df = significant_df.sort_values("Correlation", ascending=False)

    print("✅ 显著相关的目标基因：")
    print(significant_df)

    return significant_df["Gene"].tolist()


def gene_interactions(tf_target_file, significant_genes):

    # 读取 TF-Target 数据
    tf_target_df = pd.read_csv(tf_target_file, sep="\t", compression="gzip", skiprows=1)
    tf_target_df.rename(columns={"GeneSym": "source", "GeneSym.1": "target"}, inplace=True)
    # 只保留 source 和 target 都在 significant_genes 里的数据
    filtered_interactions = tf_target_df[
        (tf_target_df["source"].isin(significant_genes)) &
        (tf_target_df["target"].isin(significant_genes))
    ]

    print("✅ 目标基因之间的调控关系：")
    print(filtered_interactions[["source", "target", "weight"]])

    return filtered_interactions


def plot_network( mir100hg_targets, filtered_interactions):

    # 创建有向图
    G = nx.DiGraph()

    # **1. 添加 MIR100HG 及其 Target Genes**
    G.add_node("MIR100HG", color="red", size=1500)
    for gene in mir100hg_targets:
        G.add_node(gene, color="blue", size=800)
        G.add_edge("MIR100HG", gene)  # 连接 MIR100HG → Target Gene

    # **2. 添加 Target Genes 之间的调控关系**
    for _, row in filtered_interactions.iterrows():
        G.add_edge(row["source"], row["target"], weight=row["weight"])

    # **获取节点颜色 & 大小**
    node_colors = [G.nodes[n]["color"] for n in G.nodes]
    node_sizes = [G.nodes[n]["size"] for n in G.nodes]

    # **绘制网络图**
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=node_sizes, edge_color="gray", font_size=10, font_weight="bold")

    plt.title("MIR100HG Regulatory Network", fontsize=18)
    plt.show()



expression_file = "NTU_DATA/TOIL_Filtered_Expression_PAAD.csv"
interaction_file = "NTU_DATA/gene_attribute_edges.txt.gz"

significant_genes = filter_targets(expression_file, tf_target_df)
interactions = gene_interactions( interaction_file, significant_genes)
plot_network(significant_genes,interactions)



