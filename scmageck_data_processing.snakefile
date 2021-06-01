rule data_processing:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/scmageck_data_processing.R {config[mtx_dir]} {config[project]} {config[species]}
        """