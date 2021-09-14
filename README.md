# SCREEN
SCREEN(Single-cell CRISPR scREen data analyses and pErturbation modeliNg) is a pipeline to visualize the quality of single-cell CRISPR screens RNA-seq/ATAC-seq datasets. SCREEN has integrated three R functions: scMAGeCK_lr, Mixscape and plot function of cicero. These functions are used to estimate the regulatory score between perturbations and genes, estimate the perturbation efficiency and visualize enhancer regulatory potential, respectively.

## Dependency
	R>=4.0.3
	Seurat>=4.0.3
	scMAGeCK
	snakemake
	mixtools

## Parameter
	mtx_dir: The matrix file directory,require rds file with SeuratObject.
	sg_dir: The sgRNA file directory,include columns: cell,barcode,gene.
	prefix: The prefix directory of output file.Default current directory.
	species: The species of RNA-seq.Default Homo sapiens
	NTC: Name of Negative Control.Default NTC.
	function: Use which function in this pipeline.Default all function.
	replicate: Replicate information.Default no replicate.
	features: Features used to scale.Default high variable genes.

### Run SCREEN pipeline with one command
	snakemake --snakefile scCRISPR.snakefile --cores 4 --configfile config.yaml --config mtx_dir=input_matrix sg_dir=sgRNA_information prefix=prefix function=function_to_use
