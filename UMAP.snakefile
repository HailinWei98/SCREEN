rule UMAP:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/UMAP.R {config[rds_dir]} {config[sg_dir]}
        """