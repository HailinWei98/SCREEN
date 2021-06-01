rule data_processing:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/16_perturb_data_processing.R {config[mtx_dir]} {config[barcode_dir]} {config[gene_dir]} {config[sg_dir]} {config[project]} {config[species]}
        """