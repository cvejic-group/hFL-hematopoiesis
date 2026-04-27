#!/usr/bin/bash


# prep common variants
mkdir ~/RefData/souporcell && cd ~/RefData/souporcell
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=13aebUpEKrtjliyT9rYzRijtkNJVUk5F_" -O common_variants_grch38.vcf && rm -rf /tmp/cookies.txt
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lw4T6d7uXsm9dt39ZtEwpuB2VTY3wK1y' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lw4T6d7uXsm9dt39ZtEwpuB2VTY3wK1y" -O common_variants_hg19.vcf && rm -rf /tmp/cookies.txt
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' common_variants_grch38.vcf > common_variants_hg38.vcf

# mergebams
# merge bams of diff. libraries from the same donor
bash mergebams.sh

# run
for lib in "${pooled_libs[@]}"
do
cr=cellranger-arc202_count_${lib}_nuclei_GRCh38-2020-A-2_0_0
nohup bash souporcell2.sh -i $cr_dir/$cr/gex_possorted_bam.bam -b $cr_dir/$cr/filtered_feature_bc_matrix/barcodes.tsv.gz -k 3 -t $PPN -o $cr > ${cr}.log &
done

