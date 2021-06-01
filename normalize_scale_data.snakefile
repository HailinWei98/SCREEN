rule scale_data:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/normalize_scale_data.R {config[mtx_dir]}
        """