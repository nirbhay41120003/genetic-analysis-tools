import pandas as pd
import plotly.express as px

def plot_variants(vcf_file):
    df = pd.read_csv(vcf_file, comment="#", sep="\t", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    fig = px.scatter(df, x="POS", y="QUAL", color="FILTER", title="Variant Quality")
    fig.show()

if __name__ == "__main__":
    plot_variants("variants.vcf")
