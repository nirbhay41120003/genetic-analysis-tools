import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("gene_counts_tpm.csv")
sns.heatmap(df.iloc[:, 1:], cmap="coolwarm")

plt.title("Gene Expression Heatmap")
plt.show()
