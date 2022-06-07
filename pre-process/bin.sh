#! /bin/bash
# cat gtf/hg38_len | awk '{print $1 "\t" $NF}'
# bedtools makewindows -g gtf/hg38_len.txt -w 500 > gtf/hg38_bin.bed
# sort -k 1,1
# bedtools intersect -wa -a gtf/hg38_bin.bed -b alignment/atac.bed -c > alignment/atac.cov
#awk '$NF>0{print $1 "\t" $2 "\t" $3}' alignment/atac.cov > alignment/atac.cov.filtered
#awk '{print $4}' alignment/atac.bed > alignment/final_results/barcodes.tsv
awk 'prefix != $4 {++slab; out = sprintf($4".test", slab); prefix = $4}{print > out}' atac.bed

