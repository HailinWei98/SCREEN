rule all:
    shell:
        """
        prefix=$(pwd)
        cd {config[prefix]}
        if [ {config[function]}=='all' ];
        then
        echo "adding meta data to rds file"
        Rscript $prefix/meta_data.R {config[mtx_dir]} {config[sg_dir]} {config[species]} {config[replicate]}
        echo "finish adding meta data"
        echo "running single cell QC"
        Rscript $prefix/scQC.R {config[mtx_dir]}
        echo "finish single cell QC"
        echo "running sgRNA QC"
        Rscript $prefix/sgQC.R {config[mtx_dir]} {config[sg_dir]}
        echo "finish sgRNA QC"
        echo "running Mixscape"
        Rscript $prefix/mixscape.R {config[mtx_dir]} {config[sg_dir]} {config[NTC]}
        echo "finish Mixscape"
        echo "begin normalize and scale data"
        Rscript $prefix/normalize_scale_data.R {config[mtx_dir]} {config[features]}
        echo "finish normalize and scale"
        echo "running scMAGeCK_lr"
        Rscript $prefix/scmageck.R {config[mtx_dir]} "sg_lib.txt" {config[NTC]} {config[prefix]}
        echo "finish scMAGeCK_lr"
        echo "running UMAP"
        Rscript $prefix/UMAP.R {config[mtx_dir]} {config[sg_dir]}
        echo "all finished"
        elif [ {config[function]}=='scmageck' ];
        then
        echo "running scMAGeCK_lr"
        Rscript $prefix/scmageck.R {config[mtx_dir]} {config[sg_dir]} {config[NTC]} {config[prefix]}
        echo "finished"
        elif [ {config[function]}=='Mixscape' ];
        then
        echo "adding meta data to rds file"
        Rscript $prefix/meta_data.R {config[mtx_dir]} {config[sg_dir]} {config[species]} {config[replicate]}
        echo "finish adding meta data"
        echo "running Mixscape"
        Rscript $prefix/mixscape.R {config[mtx_dir]} {config[sg_dir]} {config[NTC]}
        echo "finished"
        elif  {config[function]}=='QC' ];
        then
        echo "running single cell QC"
        Rscript $prefix/scQC.R {config[mtx_dir]}
        echo "finish single cell QC"
        echo "running sgRNA QC"
        Rscript $prefix/sgQC.R {config[mtx_dir]} {config[sg_dir]}
        echo "finish sgRNA QC"
        elif  {config[function]}=='nc' ];
        then
        echo "begin normalize and scale data"
        Rscript $prefix/normalize_scale_data.R {config[mtx_dir]} {config[features]}
        echo "finish normalize and scale"
        fi
        """