#! /bin/bash

chr=`awk '{print $1}' alignment/atac.bed | uniq`
whole=`awk '{print $1}' gtf/hg38_len.txt | uniq`
for i in $chr
do
    if (( ${#i} > 2 )); then
        sub=`expr substr $i 1 8`
        whole_chr=`$whole | grep $sub`
        sed "s/$i/$whole_chr/g" alignment/atac.bed > alignment/atac_new.bed
    fi
done
