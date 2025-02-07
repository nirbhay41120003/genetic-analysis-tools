import numpy as np
import pandas as pd

def tpm(counts):
    rpk = counts / np.sum(counts, axis=0)
    return (rpk / np.sum(rpk, axis=0)) * 1e6

if __name__ == "__main__":
    df = pd.read_csv("gene_counts.csv")
    df_tpm = tpm(df.iloc[:, 1:])
    df_tpm.to_csv("gene_counts_tpm.csv", index=False)
