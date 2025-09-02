import pandas as pd
import numpy as np
cancer_types = ["LUAD", "PAAD", "PRAD", "SKCM", "STAD"]

# # 读取基因相互作用文件
# file_path = './NTU_DATA/gene_attribute_edges.txt.gz'
# gene_edges = pd.read_csv(file_path, sep='\t')
#
# # 筛选源基因（source）或靶基因（target）是 MIR100HG 的记录
# mir100hg_interactions = gene_edges[(gene_edges['source'] == 'MIR100HG') | (gene_edges['target'] == 'MIR100HG')]
#
# # 将筛选出的数据保存到新的 CSV 文件中
# output_file = 'NTU_DATA/gene_attribute_edge_MIR100HG.csv'
# mir100hg_interactions.to_csv(output_file, index=False)
#
# print(f"筛选出的MIR100HG相关记录已保存到 {output_file}")
#


def get_mir100hg_target_genes(edge_file):
    """
    Reads gene_attribute_edge_MIR100HG.csv to extract genes regulated by MIR100HG.
    :param edge_file: Path to the gene_attribute_edge_MIR100HG.csv file.
    :return: List of target genes affected by MIR100HG.
    """
    edge_df = pd.read_csv(edge_file)

    # 假设文件包含两列：'source'（调控基因）和 'target'（受影响基因）
    target_genes = edge_df[edge_df['source'] == 'MIR100HG']['target'].tolist()

    print(f"Found {len(target_genes)} genes regulated by MIR100HG.")
    target_genes.append("MIR100HG")
    return target_genes


def filter_expression_data(expression_file, target_genes):
    """
    Filters gene expression data to only include MIR100HG and its target genes.
    :param expression_file: Path to the TOIL_RSEM_TPM_XXXX.csv file.
    :param target_genes: List of genes to retain.
    :return: Filtered DataFrame.
    """
    expression_df = pd.read_csv(expression_file, index_col="HGNC_symbol")

    # 仅保留目标基因
    filtered_df = expression_df.loc[expression_df.index.intersection(target_genes)]

    print(f"Filtered gene expression data shape: {filtered_df.shape}")
    return filtered_df

def save_filtered_expression_data(df, save_path):
    """
    Saves the filtered expression data to a CSV file.
    :param df: Filtered DataFrame.
    :param save_path: Path to save the file.
    """
    df.to_csv(save_path)
    print(f"Filtered expression data saved to {save_path}")


# # 读取 MIR100HG 影响的基因
# edge_file = "./NTU_DATA/gene_attribute_edge_MIR100HG.csv"
# mir100hg_targets = get_mir100hg_target_genes(edge_file)
#
# print("Target gene list:", mir100hg_targets)
#

# -------------------------------------------


# 示例：处理 LUAD 数据
# expression_file = "./NTU_DATA/TOIL_RSEM_TPM_STAD.csv"
# filtered_expression_df = filter_expression_data(expression_file, mir100hg_targets)
#
# # 预览数据
# print(filtered_expression_df.head())
#
#
# # 保存数据
# save_path = "./NTU_DATA/TOIL_Filtered_Expression_STAD.csv"
# save_filtered_expression_data(filtered_expression_df, save_path)


# -------------------------------------------
"""寻找探针对应的CpG 点位"""

def get_mir100hg_regions(annotation_file):
    """
    Extracts all genomic regions of MIR100HG from the gene annotation file.
    :param annotation_file: Path to the gene annotation file.
    :return: DataFrame with chromosome, start, and end positions of MIR100HG.
    """
    # 读取基因注释文件
    annotation_df = pd.read_csv(annotation_file, delimiter="\t")

    # 筛选所有 MIR100HG 相关的行
    mir100hg_regions_df = annotation_df[(annotation_df['symbol'] == 'MIR100HG' )& (annotation_df['type'] == 'hg19_genes_promoters')]

    print(f"Found {mir100hg_regions_df.shape[0]} regions for MIR100HG.")
    print(mir100hg_regions_df)
    return mir100hg_regions_df


def filter_mir100hg_probes(probe_map_file, mir100hg_regions_df):
    """
    Filters CpG probes located within all MIR100HG genomic regions.
    :param probe_map_file: Path to the probe mapping file.
    :param mir100hg_regions_df: DataFrame with MIR100HG genomic regions.
    :return: DataFrame with filtered probes.
    """
    # 读取 probeMap 文件
    probe_map_df = pd.read_csv(probe_map_file, delimiter="\t")

    # 存储所有符合条件的探针
    all_mir100hg_probes = pd.DataFrame()

    # 遍历每个 MIR100HG 相关的基因区间
    for _, row in mir100hg_regions_df.iterrows():
        chrom, start, end = row['chr'], row['start'], row['end']

        # 筛选出匹配该基因区间的 CpG 探针
        filtered_probes = probe_map_df[
            (probe_map_df['gene'] == "MIR100HG") &
            (probe_map_df['chrom'] == chrom) &  # 染色体匹配
            (probe_map_df['chromStart'] >= start) &  # 起始位置在范围内
            (probe_map_df['chromEnd'] <= end)  # 终止位置在范围内
        ]

        # 合并结果
        all_mir100hg_probes = pd.concat([all_mir100hg_probes, filtered_probes])

    # 去重（避免多个区间重复探针）
    all_mir100hg_probes = all_mir100hg_probes.drop_duplicates()

    print(f"Found {all_mir100hg_probes.shape[0]} unique CpG probes in MIR100HG regions.")
    print(all_mir100hg_probes)
    return all_mir100hg_probes

# # 读取 MIR100HG 相关区间，并筛选 MIR100HG 相关的探针
# annotation_file = "./NTU_DATA/geneAnnotation_hg19_basicgenes.txt"
# probe_map_file = "./NTU_DATA/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy"
#
# mir100hg_regions_df = get_mir100hg_regions(annotation_file)
# mir100hg_probes_df = filter_mir100hg_probes(probe_map_file, mir100hg_regions_df)
#
# mir100hg_probes_df.to_csv("./NTU_DATA/probeMap_MIR100HG_related.csv", sep= "\t")
#

# -------------------------------------------
"""进一步筛选五种癌症中的甲基化点位"""

def filter_mir100hg_methylation(mir100hg_probes_file, mir100hg_methylation_file):
    """
    Filters the MIR100HG-related methylation data, keeping only probes listed in probeMap_MIR100HG_related.csv.
    :param mir100hg_probes_file: Path to the probe map file (only relevant probes).
    :param mir100hg_methylation_file: Path to the methylation data file to be filtered.
    """

    # 读取 probeMap 文件
    mir100hg_probe_map = pd.read_csv(mir100hg_probes_file, delimiter='\t')

    # 获取 MIR100HG 相关的探针
    mir100hg_probes = mir100hg_probe_map['#id'].tolist()

    # 读取 MIR100HG_Methylation_XXXX.csv
    methylation_df = pd.read_csv(mir100hg_methylation_file, delimiter=",", index_col=0)

    # 过滤甲基化数据（匹配索引，而不是列）
    mir100hg_methylation_df = methylation_df.loc[methylation_df.index.intersection(mir100hg_probes)]

    # 保存并覆盖原文件
    mir100hg_methylation_df.to_csv(mir100hg_methylation_file, sep="\t", index=False)
    print(f"Updated {mir100hg_methylation_file}: {mir100hg_methylation_df.shape[0]} probes retained.")

# # 进一步处理 Methylation数据
# for cancer in cancer_types:
#     methylation_file = f"NTU_DATA/MIR100HG_Methylation_{cancer}.csv"
#     filter_mir100hg_methylation("NTU_DATA/probeMap_MIR100HG_related.csv", methylation_file)

# -------------------------------------------
"""----------------------------------data-fusion-----------------------------"""

# 数据格式转化 TCGA-XX-XXXX => TCGA_XX_XXXX
# def format_sample_id(sample_id):
#     """
#     Converts Methylation sample ID (TCGA-XX-XXXX-YY) to Expression format (TCGA_XX_XXXX).
#     """
#     return "_".join(sample_id.split("-")[:3])  # 只保留前三段
#
# for cancer in cancer_types:
#     # 读取 LUAD 甲基化数据
#     data_folder = "./NTU_DATA"
#     methylation_df = pd.read_csv(f"./NTU_DATA/MIR100HG_Methylation_{cancer}.csv", index_col=0)
#     methylation_df.columns = [format_sample_id(sample) for sample in methylation_df.columns]
#     methylation_df.to_csv(f"./NTU_DATA/MIR100HG_Methylation_{cancer}.csv", sep="\t", index=True)
#

def merge_mir100hg_expression_to_methylation(methylation_file, expression_file, output_file):
    """
    Adds a new row to the methylation file containing MIR100HG expression values.
    :param methylation_file: Path to the methylation data file.
    :param expression_file: Path to the gene expression data file.
    :param output_file: Path to save the updated methylation file.
    """
    # **读取 Methylation 甲基化数据**
    methylation_df = pd.read_csv(methylation_file, sep="\t", index_col=0)

    # **读取 Expression 基因表达数据**
    expression_df = pd.read_csv(expression_file, index_col="HGNC_symbol")

    mir100hg_expression = expression_df.loc["MIR100HG"]

    # **匹配样本**
    common_samples = methylation_df.columns.intersection(mir100hg_expression.index)

    print(f"✅ Total Methylation samples: {len(methylation_df.columns)}")
    print(f"✅ Total Expression samples: {len(expression_df.columns)}")
    print(f"✅ Matching samples found: {len(common_samples)}")

    # **创建 MIR100HG 新行**
    mir100hg_row = pd.Series(np.nan, index=methylation_df.columns, name="MIR100HG")
    mir100hg_row[common_samples] = mir100hg_expression[common_samples]

    # **添加到 Methylation 数据**

    methylation_df = pd.concat([methylation_df, mir100hg_row.to_frame().T])

    # **去除包含 NaN（空值）的列**
    methylation_df.dropna(axis=1, inplace=True)

    # **保存数据**
    methylation_df.to_csv(output_file, sep=",",index=True)
    print(f"✅ Updated file saved: {output_file}")
#
# for cancer in cancer_types:
#     merge_mir100hg_expression_to_methylation(
#         f"./NTU_DATA/MIR100HG_Methylation_{cancer}.csv",
#         f"./NTU_DATA/TOIL_Filtered_Expression_{cancer}.csv",
#         f"./NTU_DATA/MIR100HG_Methylation_{cancer}.csv"
#     )
#
