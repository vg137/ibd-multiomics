# Integrating GWAS and Gene-Set Expression in IBD

This repository contains a three-notebook analysis of inflammatory bowel disease
(IBD). It combines genome-wide association study (GWAS) summary statistics with
bulk RNA-seq expression data and Hallmark gene sets to identify IBD-associated
genes and rank biological programs that contain them.

The expression dataset is bulk RNA-seq, not single-cell RNA-seq. The notebooks
include analysis code and narrative, but the input data are not included in
this repository. The helper functions used by Notebook 1 are provided in
[`nb1_functions.py`](nb1_functions.py).

## Workflow

Run the notebooks in this order:

1. [`nb1_snp2genes.ipynb`](nb1_snp2genes.ipynb) downloads or loads the IBD GWAS
	summary statistics, performs basic quality-control checks, filters SNPs at
	$p < 5 \times 10^{-8}$, and assigns each significant SNP to its nearest gene
	transcription start site using GENCODE annotations.
2. [`nb2_gene_set_expressions.ipynb`](nb2_gene_set_expressions.ipynb) loads
	RNA-seq TPM values from GEO series GSE57945, maps NCBI gene IDs to symbols,
	and calculates mean expression for each MSigDB Hallmark gene set in each
	sample.
3. [`nb3_integrating_gwas_and_expressions.ipynb`](nb3_integrating_gwas_and_expressions.ipynb)
	intersects the GWAS-derived genes with the expressed Hallmark gene sets,
	ranks gene sets by the size of that intersection, and identifies genes that
	occur in multiple relevant gene sets.

## Data sources

| Resource | Dataset or version | Used for |
| --- | --- | --- |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | GCST90292538 | IBD summary statistics |
| [GENCODE](https://www.gencodegenes.org/human/) | GRCh38.p14 | Gene coordinates and transcription start sites |
| [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57945) | GSE57945 | Bulk RNA-seq TPM expression values |
| [NCBI Gene](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) | `Homo_sapiens.gene_info.gz` | Gene ID to symbol mapping |
| [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) | Hallmark, `h.all.v2025.1.Hs.symbols.gmt` | Curated gene sets |

The exact files must be obtained from their respective providers. Check each
provider's terms and citation requirements before downloading or redistributing
data.

## Requirements

- Python 3 with Jupyter or Google Colab
- `pandas`, `numpy`, `matplotlib`, `seaborn`, `pyarrow`, and `nbformat`
- The `nb1_functions` module imported by Notebook 1
- Enough local or Google Drive storage for the GWAS and annotation parquet files

The notebooks were written for Google Colab and mount Google Drive at the start
of each workflow. When running locally, keep `nb1_functions.py` in the same
directory as Notebook 1 so that `import nb1_functions` resolves. In Colab,
clone or upload this repository and add its directory to the Python import path
before running Notebook 1. The repository is not currently packaged as a
standalone Python environment and does not contain `requirements.txt`.

## Expected data layout

Set the notebook paths to a project directory containing this structure. Raw
inputs are provider downloads; parquet files and analysis outputs are generated
by the notebooks.

```text
ibd_data/
├── gwas/raw/GCST90292538_harmonised/
│   ├── GCST90292538.h.tsv.gz
│   └── GCST90292538.h.tsv.gz-meta.yaml
├── gwas/parquet/GCST90292538.parquet
├── annotations/raw/gencode.v49.annotation.gtf.gz
├── annotations/parquet/gencode.v49.annotation.parquet
├── expressions/GSE57945/GSE57945_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz
├── geneID2symbol/Homo_sapiens.gene_info.gz
├── gene_sets/h.all.v2025.1.Hs.symbols.gmt
├── snp2gene/GCST90292538_snp2gene.parquet
├── gene_set_expressions/hallmark_expressions.parquet
├── gwas_integrated_expressions/ranked_intersection.parquet
└── gene_support_in_gene_sets/gene_support.parquet
```

Notebook 1 creates intermediate GWAS and annotation parquet files and the
SNP-to-gene mapping. Notebook 2 creates the Hallmark gene-set expression file.
Notebook 3 consumes those two files and writes the integration outputs shown
above. Commented-out save steps in the notebooks must be enabled when a file
does not yet exist.

## Results and limitations

The notebooks are exploratory analysis workflows, not a validated causal-
inference pipeline. A nearest-gene assignment does not establish that a gene is
causal, and ranking a gene set by overlap does not test statistical enrichment.
The current notebooks report candidate genes and inflammation- or immune-related
Hallmark programs; reproduce the notebooks and inspect their plots and saved
tables before treating those observations as final results.
