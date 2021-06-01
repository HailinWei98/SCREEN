rule scmageck:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        Rscript $prefix/scmageck.R {config[mtx_dir]} {config[sg_dir]} {config[NTC]} {config[prefix]}
        """