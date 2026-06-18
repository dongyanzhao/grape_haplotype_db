# 2025.10.20
# Some trait markers for grape dartag are not on PN40024 reference genome
# These trait markers are mostly from other grape varieties
# There are lots of IUPAC nucleotide codes and indels (-) in the marker sequences
# Some trait markers have "Other" type of microhaplotypes prevalent in the validation results
# Therefore, they need to be processed separately from the chromosome-based markers

# Extract MADC for these markers
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db02_update_snpID_in_madc_v1.py \
/Users/dz359/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_nonChr_lut.csv \
/Users/dz359/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10633_Lance_4plates/OrderAppendix_1_DG25-10633_flipRefAlt_nonChr.csv
```


# ========================================
# Use multiple sequence alignment to refine the ref/alt allele sequences for these markers
# ========================================
# 1. All alleles alignment
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_madc_multi_sequence_aln.py \
/Users/dz359/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10633_Lance_4plates/OrderAppendix_1_DG25-10633_flipRefAlt_nonChr_snpID.csv \
/Users/dz359/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/nonChr_marker_aln \
/Users/dz359/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20251021-nonChr_allAlleles \
-first_bp=125
```

# 2. Ref vs. Alt only alignment
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_madc_multi_sequence_aln.py \
/Users/dz359/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10633_Lance_4plates/OrderAppendix_1_DG25-10633_flipRefAlt_nonChr_snpID.csv \
/Users/dz359/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/nonChr_marker_aln \
/Users/dz359/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20251021-nonChr_refAlt \
-ref_alt_only=Y \
-aln_suffix=refAlt \
-first_bp=125
```

# 3. Get Ref/Alt sequences based on sequence alignments above
## This is done manually.
## /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20251021-nonChr_ref_alt_rev.fa

