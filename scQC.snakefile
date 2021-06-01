rule scQC:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/scQC.R {config[mtx_dir]}
        """