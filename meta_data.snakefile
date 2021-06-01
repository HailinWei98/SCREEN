rule meta_data:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/meta_data.R {config[mtx_dir]} {config[sg_dir]} {config[species]}
        """