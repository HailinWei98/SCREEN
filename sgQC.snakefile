rule sgQC:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/sgQC.R {config[mtx_dir]} {config[sg_dir]} {config[project]}
        """