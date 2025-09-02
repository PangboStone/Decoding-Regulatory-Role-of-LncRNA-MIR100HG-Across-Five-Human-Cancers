import pandas as pd
import os
cancer_types = ["LUAD", "PAAD", "PRAD", "SKCM", "STAD"]

def read_methylation_data(data_folder, cancer_type="LUAD"):
    """
    Reads the DNA methylation data for the specified cancer type.
    :param data_folder: Path to the folder containing the data files.
    :param cancer_type: The specific cancer type (e.g., LUAD, PAAD, etc.).
    :return: DataFrame containing methylation data.
    """
    file_name = f"TCGA.{cancer_type}.sampleMap_HumanMethylation450.gz"
    file_path = os.path.join(data_folder, file_name)

    if os.path.exists(file_path):
        df = pd.read_csv(file_path, compression='gzip', delimiter='\t')
        return df
    else:
        print(f"Warning: {file_name} not found.")
        return None


def filter_mir100hg_methylation(methylation_df, probe_map_file):
    """
    Filters the methylation dataset to only include MIR100HG-related CpG probes.
    :param methylation_df: DataFrame with all CpG probe methylation data.
    :param probe_map_file: Path to the probe mapping file.
    :return: Filtered DataFrame with only MIR100HG-related CpG probes.
    """
    # 读取 probeMap 文件
    probe_map = pd.read_csv(probe_map_file, delimiter='\t')

    # 获取 MIR100HG 相关的探针
    mir100hg_probes = probe_map[probe_map['gene'] == 'MIR100HG']['#id'].tolist()

    # 过滤甲基化数据
    mir100hg_methylation_df = methylation_df[methylation_df['sample'].isin(mir100hg_probes)]

    return mir100hg_methylation_df


def save_methylationData_to_csv(df, save_path):
    """
    Saves the filtered methylation data to a CSV file.
    :param df: DataFrame to save.
    :param save_path: Path to save the file.
    """
    df.to_csv(save_path, index=False)
    print(f"Filtered methylation data saved to {save_path}")

# 进一步处理 Methylation数据
for cancer in cancer_types:
    # 读取 LUAD 甲基化数据
    data_folder = "./NTU_DATA"
    probe_map_file = "./NTU_DATA/probeMap_MIR100HG_related.csv"
    methylation_df = read_methylation_data(data_folder, cancer)
    mir100hg_methylation_df = filter_mir100hg_methylation(methylation_df, probe_map_file)
    save_methylationData_to_csv(mir100hg_methylation_df, f"./NTU_DATA/MIR100HG_Methylation_{cancer}.csv")

    # print(f"Filtered methylation data shape: {mir100hg_methylation_df.shape}")
    # print(mir100hg_methylation_df.head())

