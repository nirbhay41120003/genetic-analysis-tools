import plotly.express as px
import pandas as pd

df = pd.read_csv("variants.vcf", comment="#", sep="\t", names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
fig = px.scatter(df, x="POS", y="QUAL", title="Variant Quality", labels={"POS": "Genomic Position", "QUAL": "Quality Score"})
fig.show()
