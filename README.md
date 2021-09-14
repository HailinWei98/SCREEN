# SCREEN
SCREEN(<strong>S</strong>ingle-cell <strong>C</strong>RISPR sc<strong>RE</strong>en data analyses and p<strong>E</strong>rturbation modeli<strong>N</strong>g) is a pipeline to visualize the quality of single-cell CRISPR screens RNA-seq/ATAC-seq datasets. SCREEN has integrated three R functions: scMAGeCK_lr, Mixscape and plot function of cicero. These functions are used to estimate the regulatory score between perturbations and genes, estimate the perturbation efficiency and visualize enhancer regulatory potential, respectively.

## Dependency
	R >= 4.0.3
	Seurat >= 4.0.3
	scMAGeCK
	snakemake
	mixtools

## Parameter
<center><code>mtx_dir</code> : The matrix file directory, require rds file with SeuratObject.</center><br>
<center><code>sg_dir</code> : The sgRNA file directory, include columns: cell, barcode, gene.</center><br>
<center><code>prefix</code> : The prefix directory of output file. Default current directory.</center><br>
<center><code>species</code> : The species of RNA-seq. Default Homo sapiens.</center><br>
<center><code>NTC</code> : Name of Negative Control. Default 'NTC'.</center><br>
<center><code>function</code> : Use which function in this pipeline. Default all function.</cencer><br>
<center><code>replicate</code> : Replicate information. Default no replicate.</cencer><br>
<center><code>features</code> : Features used to scale. Default high variable genes.</center>

### Run SCREEN pipeline with one command
	snakemake --snakefile scCRISPR.snakefile --cores 4 --configfile config.yaml --config mtx_dir=input_matrix sg_dir=sgRNA_information prefix=prefix function=function_to_use
