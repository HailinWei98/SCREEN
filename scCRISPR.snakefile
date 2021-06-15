rule all:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/meta_data.R {config[mtx_dir]} {config[sg_dir]} {config[species]}
        Rscript $prefix/scQC.R {config[mtx_dir]}
        Rscript $prefix/sgQC.R {config[mtx_dir]} {config[sg_dir]} {config[project]}
        Rscript $prefix/mixscape.R {config[mtx_dir]} {config[sg_dir]}
        Rscript $prefix/normalize_scale_data.R {config[mtx_dir]}
        Rscript $prefix/scmageck.R {config[mtx_dir]} {config[sg_dir]} {config[NTC]} {config[prefix]}
        Rscript $prefix/UMAP.R {config[rds_dir]} {config[sg_dir]}
        """