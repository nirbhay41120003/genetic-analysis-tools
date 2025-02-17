import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_expression(data_file):
    df = pd.read_csv(data_file)
    plt.figure(figsize=(10, 5))
    sns.heatmap(df.iloc[:, 1:], cmap="coolwarm")
    plt.title("Gene Expression Heatmap")
    plt.show()

if __name__ == "__main__":
    plot_expression("gene_expression/gene_counts_tpm.csv")
