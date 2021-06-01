rule data_processing:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/18_crop_data_processing.R {config[mtx_dir]} {config[sg_dir]} {config[project]} {config[species]}
        """