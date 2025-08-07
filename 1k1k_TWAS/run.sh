mkdir -p data/intermediate/sumstats
tail -n+2 data/gwas/gwas_metadata.txt | cut -f2 | head -3 > data/intermediate/traits.3.txt
cat data/intermediate/traits.3.txt | while read trait; do
    echo $trait
    # python ~/tools/ldsc_py3/munge_sumstats.py \
    #     --sumstats data/gwas/imputed_gwas_hg38_1.1/imputed_${trait}.txt.gz \
    #     --out data/intermediate/sumstats/${trait}
    echo -e "SNP\tA1\tA2\tZ" > data/intermediate/sumstats/${trait}.txt
    zcat data/gwas/imputed_gwas_hg38_1.1/imputed_${trait}.txt.gz \
        | tail -n+2 \
        | cut -f2,5,6,10 \
        >> data/intermediate/sumstats/${trait}.txt
done

## Convert IDs in LDREF to match genotypes used by Pantry
for chrom in {1..22}; do
    echo $chrom
    awk 'BEGIN {OFS="\t"} NR==FNR {id[$1]=$2; next} {if ($2 in id) $2=id[$2]; else $2="NA"; print}' \
        <(zcat data/gwas/imputed_gwas_hg38_1.1/imputed_UKB_1160_Sleep_duration.txt.gz) \
        data/LDREF/1000G.EUR.${chrom}.bim \
        > data/LDREF_b38ids/1000G.EUR.${chrom}.bim
    mv -i data/LDREF/1000G.EUR.${chrom}.bed data/LDREF_b38ids/
    mv -i data/LDREF/1000G.EUR.${chrom}.fam data/LDREF_b38ids/
done

## Subset Pantry genotypes to SNPs in LDREF
cat data/LDREF_b38ids/1000G.EUR.*.bim | cut -f2 | sort | uniq > data/intermediate/LD_SNPs.txt
