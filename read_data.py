import os
import pandas as pd

def read_gene_expression_data(data_folder):
    """
    Reads gene expression data from the NTU_DATA folder.
    :param data_folder: Path to the folder containing the data files.
    :return: Dictionary with cancer types as keys and DataFrames as values.
    """
    gene_expression_files = [
        "TOIL RSEM TPM LUAD.CSV",
        "TOIL RSEM TPM PAAD.CSV",
        "TOIL RSEM TPM PRAD.CSV",
        "TOIL RSEM TPM SKCM.CSV",
        "TOIL RSEM TPM STAD.CSV"
    ]

    gene_expression_data = {}

    for file_name in gene_expression_files:
        file_path = os.path.join(data_folder, file_name)
        if os.path.exists(file_path):
            cancer_type = file_name.split(".")[2]
            df = pd.read_csv(file_path)
            gene_expression_data[cancer_type] = df
        else:
            print(f"Warning: {file_name} not found in the directory.")

    return gene_expression_data


def read_methylation_data(data_folder):
    """
    Reads methylation data from the NTU_DATA folder.
    :param data_folder: Path to the folder containing the data files.
    :return: Dictionary with cancer types as keys and DataFrames as values.
    """
    methylation_files = [
        "TCGA.LUAD.sampleMap HumanMethylation450.gz",
        "TCGA.PAAD.sampleMap HumanMethylation450.gz",
        "TCGA.PRAD.sampleMap HumanMethylation450.gz",
        "TCGA.SKCM.sampleMap HumanMethylation450.gz",
        "TCGA.STAD.sampleMap HumanMethylation450.gz"
    ]

    methylation_data = {}

    for file_name in methylation_files:
        file_path = os.path.join(data_folder, file_name)
        if os.path.exists(file_path):
            cancer_type = file_name.split(".")[1]
            df = pd.read_csv(file_path, compression='gzip')
            methylation_data[cancer_type] = df
        else:
            print(f"Warning: {file_name} not found in the directory.")

    return methylation_data


def read_transcription_factor_data(data_folder):
    """
    Reads transcription factor and target gene association data.
    :param data_folder: Path to the folder containing the data files.
    :return: DataFrame containing the transcription factor data.
    """
    file_path = os.path.join(data_folder, "gene attribute edges.txt.gz")
    if os.path.exists(file_path):
        return pd.read_csv(file_path, compression='gzip', delimiter='\t', header=None)
    else:
        print("Warning: Transcription factor data file not found.")
        return None


def read_clinical_data(data_folder):
    """
    Reads clinical data from the NTU_DATA folder.
    :param data_folder: Path to the folder containing the data files.
    :return: DataFrame containing the clinical information.
    """
    file_path = os.path.join(data_folder, "Survival SupplementalTable S1 20171025 xena sp")
    if os.path.exists(file_path):
        return pd.read_csv(file_path)
    else:
        print("Warning: Clinical data file not found.")
        return None


def read_all_data(data_folder="./NTU_DATA"):
    """
    Reads all necessary data from the NTU_DATA folder.
    :param data_folder: Path to the folder containing the data files.
    :return: Dictionary with all datasets loaded.
    """
    data = {
        "gene_expression": read_gene_expression_data(data_folder),
        "methylation": read_methylation_data(data_folder),
        "transcription_factors": read_transcription_factor_data(data_folder),
        "clinical_data": read_clinical_data(data_folder)
    }

    return data

if __name__ == "__main__":
    # Test the data loading functionality
    data = read_all_data()
    print("Data loading complete.")
    for key, value in data.items():
        if isinstance(value, dict):
            print(f"{key}: {len(value)} datasets loaded.")
        elif value is not None:
            print(f"{key}: {value.shape[0]} rows loaded.")


"""
对于数据probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy：
#id：探针的唯一标识符，例如 cg13332474。
gene：与探针相关的基因符号。如果有多个基因，它们会用逗号分隔。例如，cg00651829 与 RSPH14 和 GNAZ 相关。
chrom：探针所在的染色体，例如 chr7 表示第7号染色体。
chromStart 和 chromEnd：探针在染色体上的起始和结束位置，表示探针覆盖的基因组区域。例如，cg13332474 的起始位置是 25935146，结束位置是 25935148。
strand：表示探针所在的DNA链，通常是正链（+）或负链（-）。在你的数据中，这一列可能为空（用 . 表示）。
"""
