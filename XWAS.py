
from pyspark.sql import SparkSession
import hail as hl
import os

builder = (
    SparkSession
    .builder
    .enableHiveSupport()
)
spark = builder.getOrCreate()
hl.init(sc=spark.sparkContext)
###.idx

bgen_path = "/mnt/project/Bulk/Imputation/UKB imputation from genotype"


filename='ukb22828_cX_b0_v3.bgen'
file_url = f"file://{bgen_path}/{filename}"
hl.index_bgen(path=file_url,
              index_file_map={file_url:f"hdfs:///{filename}.idx2"},
              reference_genome="GRCh37",
              contig_recoding={"0X": "X"},
              skip_invalid_loci=False)



index_file_map = {}
index_file_map[f"file://{bgen_path}/{filename}"] = f"hdfs:///{filename}.idx2"
##输入mt
mt = hl.import_bgen(path=f"file://{bgen_path}/ukb22828_cX_b0_v3.bgen", # regex can be used if genomic data is in multiple BGEN files
                    entry_fields=['GT', 'GP'],
                    sample_file=f"file://{bgen_path}/ukb22828_cX_b0_v3.sample",
                    n_partitions=None,
                    block_size=None,
                    index_file_map=index_file_map,
                    variants=None,)

table = (hl.import_table('file:///mnt/project/Multiple_myeloma/pheno_cov.txt', impute=True, types={'eid': hl.tstr})
         .key_by('eid'))

mt = mt.annotate_cols(pheno=table[mt.s])

#female--------------------------------------------
gwas = hl.linear_regression_rows(
    y=mt.pheno.female_Multiple_myeloma,
    x=mt.GT.n_alt_alleles(),
    covariates=[0.0,mt.pheno.age,mt.pheno.pca_1,mt.pheno.pca_2,mt.pheno.pca_3,mt.pheno.pca_4,mt.pheno.pca_5,mt.pheno.pca_6,mt.pheno.pca_7,
mt.pheno.pca_8,mt.pheno.pca_9,mt.pheno.pca_10,mt.pheno.pca_11,mt.pheno.pca_12,mt.pheno.pca_13,mt.pheno.pca_14,mt.pheno.pca_15,mt.pheno.pca_16,mt.pheno.pca_17,mt.pheno.pca_18,
mt.pheno.pca_19,mt.pheno.pca_20,mt.pheno.age2]
)

gwas2=gwas.row
gwas2.export('female_Multiple_myeloma_chrX.tsv.gz', delimiter='\t')

%%bash
hdfs dfs -get female_Multiple_myeloma_chrX.tsv.gz
dx upload female_Multiple_myeloma_chrX.tsv.gz


#male--------------------------------------------
gwas = hl.linear_regression_rows(
    y=mt.pheno.male_Multiple_myeloma,
    x=mt.GT.n_alt_alleles(),
    covariates=[0.0,mt.pheno.age,mt.pheno.pca_1,mt.pheno.pca_2,mt.pheno.pca_3,mt.pheno.pca_4,mt.pheno.pca_5,mt.pheno.pca_6,mt.pheno.pca_7,
mt.pheno.pca_8,mt.pheno.pca_9,mt.pheno.pca_10,mt.pheno.pca_11,mt.pheno.pca_12,mt.pheno.pca_13,mt.pheno.pca_14,mt.pheno.pca_15,mt.pheno.pca_16,mt.pheno.pca_17,mt.pheno.pca_18,
mt.pheno.pca_19,mt.pheno.pca_20,mt.pheno.age2]
)

gwas2=gwas.row
gwas2.export('male_Multiple_myeloma_chrX.tsv.gz', delimiter='\t')

%%bash
hdfs dfs -get male_Multiple_myeloma_chrX.tsv.gz
dx upload male_Multiple_myeloma_chrX.tsv.gz

