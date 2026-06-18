# 2025.9.4
# Building the grape reference and alternative allele database
# ===========! IMPORTANT ===========
# 1. Markers with the same DArTag amplicons
# There are four pairs of duplicate sequences
AGL11Exon7_000000071	AGL11Exon7_000000109
chr13_016075743	chr13_016075765
chr13_027057789	chr13_027057798
chr17_017142764	chr17_017142822
# SOLUTION:
## a. Remove the second markers in the duplicate pairs from marker LUT and the final RefAlt database

# 2. Markers from Muscadine were prepared using Stacks 
# The `populations` program reads the variant calls and outputs a VCF file. While the initial variant detection is reference-based, the final REF allele reported in the VCF may not always correspond directly to the reference genome's base.For a given locus, the reference allele in the VCF output is determined by the allele with the highest count in the first individual processed by cstacks.
# This behavior of using the most common allele as the REF can cause the REF allele to be "flipped" relative to the actual reference genome, especially if a derived allele is common in your population.
# SOLUTION:
## a. I will correct the Ref, Alt, RefMatch, and AltMatch IDs in raw MADC file based on the actual reference genome sequence.
## b. Prepare Ref and Alt flanking sequences based on the reference genome
## c. Prepare a lookup table (LUT) to link the original marker IDs to the updated Ref and Alt alleles and their flanking sequences.
## d. Provide this information to researchers needing it

# 3. Markers from Muscadine are inconsistent between probe design and the actual flanking sequences extracted from the reference genome.
# Checked the Ref/Alt matches those from the original Stacks vcf
# Stacks report Ref/Alt with majority rule, which may not always be the same as the actual REF/ALT reported in the VCF.

# SOLUTION:
## a. I will remove these markers from the probe design file and the final RefAlt database
non-matching marker:	Ref in flank seq	[ref, alt] in stacks vcf
chr15_020072353	T	['G', 'A']
chr18_029482967	A	['G', 'C']
chr18_034467440	A	['T', 'C']


# 4. Markers not reported in DArTag report
chr12_017277339
chr13_023570339
chr15_000741434
chr16_010525391
chr19_001417775
# SOLUTION:
## a. I will remove these markers from the probe design file and the final RefAlt database

# 5. Markers not from nonRef sequences contain hyphens in the flanking sequences and they ("-") will be problematic for downstream analysis. Additionally, the RefAlt defined in the probe design file may not match the actual bases at the target SNP position in the flanking sequences. I suspect these issues are due to the use of random allele sequences selected among several alleles identified/collected by Cheng
# SOLUTION:
## a. I will remove hyphen.
## b. Adjust the target marker positions in the lut file accordingly.
## c. Prepare Ref and Alt flanking sequences based on the reference genome or non-reference sequences separately.

# 5. Trait and QTL markers with no similar sequences in the reference genome
## P1_apt3_7641402_1002_VariantMasked_000000112
## TA_1_7644532_1002_VariantMasked_000000121

# ===========! IMPORTANT ===========


# ======================
# Detailed analyses and steps
# ======================
# Check if there are duplicate sequences
```bash
py /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/02_fixedAlleleID/step00_check_allele_uniqueness_AND_update_madc_v1.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10633_Vvinifera/OrderAppendix_1_DG25-10633.csv \
17
```
# There are four pairs of duplicate sequences
AGL11Exon7_000000071	AGL11Exon7_000000109
chr13_016075743	chr13_016075765
chr13_027057789	chr13_027057798
chr17_017142764	chr17_017142822
>> Remove the 4 duplicated loci manually

# ===================
# Processing reference-based markers
# ===================
# 1. Prepare LUT from probe design file manually 
# <20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut.csv>

# 2. Prepare sfetch keys based on snpID lut
```bash
py /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db01_get_reference_sfetch_keys_from_snpID_lut_v1.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut.csv \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/ref/PN40024/PN40024_assembly_v2_upper_rename_nonRef.fa.len \
300
```
# Total records written out:  2984


# 3. Index reference genome
```bash
esl-sfetch --index /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/ref/PN40024/PN40024_assembly_v2_upper_rename_nonRef.fa
```


# 4. Extract flanking sequences
```bash
esl-sfetch -Cf /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/ref/PN40024/PN40024_assembly_v2_upper_rename_nonRef.fa \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys.txt \
> /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys.fa
```


# ==============
# Markers from Muscadine are inconsistent between probe design and the actual flanking sequences extracted from the reference genome.
# Checked the Ref/Alt matches those from the original Stacks vcf
# Stacks report Ref/Alt with majority rule, which may not always be the same as the actual REF/ALT reported in the VCF.
# ==============
## Check Ref and Alt amplicon sequences from MADC report
### Extract Ref and Alt sequences from MADC report
```bash
python3 /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db03_ext_ref_alt_amp_AND_gen_match_lut_from_madc_v1.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10632_Margaret_6plates/DG25-10632_MADC.csv
```

### Conduct BLAST alignment
```bash
blastn -task blastn-short -dust no -soft_masking false -word_size 7 \
-db /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/ref/PN40024/PN40024_assembly_v2_upper_rename.fa \
-query /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/grape_allele_db_v001.fa \
-out /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/grape_allele_db_v001.fa_PN40024.bn \
-evalue 1e-5 -num_threads 6 -max_target_seqs 200 -max_hsps 5 \
-outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
```

### Check if Alt allele matches 100% to the reference genome
```bash
awk -F'\t' '$1 ~ /Alt_0002/ && $9==109 && $10 == 100 && $11 == 100' \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10632_Margaret_6plates/DG25-10632_MADC_ref_alt_amplicons.fa_PN40024.bn > /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10632_Margaret_6plates/DG25-10632_MADC_ref_alt_amplicons.fa_PN40024.bn.alt 
```



## Prepare 601 bp Ref and Alt flanking sequences based on RefAlt in LUT
```bash
py /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db01_prep_ref_alt_flankSeq_from_lut_v1.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut.csv \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys.fa \
300
```
[INFO] Number of markers with Ref and Alt swapped:  119


## Add passport to the markers with Ref and Alt swapped
```bash
python3 /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/04_util/util_add_passport_to_pca.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys_flippedRefAlt.csv \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20240812-Grape_3K_DArTag_probe.csv \
--index_colID MarkerID
```
## All markers are from Muscadine



## Confirm the 601 bp Ref/Alt sequences by BLAST the sequences to the reference genome PN40024
```bash
blastn -task blastn-short -dust no -soft_masking false -word_size 7 \
-db /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/ref/PN40024/PN40024_assembly_v2_upper_rename.fa \
-query /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys_ref_alt.fa \
-out /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys_ref_alt.fa_PN40024.bn \
-evalue 1e-5 -num_threads 6 -max_target_seqs 5 \
-outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
```


# make blast database
```bash
makeblastdb -in /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys_ref_alt.fa -dbtype nucl
```


# 6. Prepare MADC file with updated Ref, Alt, RefMatch, and AltMatch
```bash
py /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db01_flip_refAlt_in_raw_MADC.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys_flippedRefAlt.csv \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10632_Margaret_6plates/DG25-10632_MADC.csv
```
# Updated 119 markers:


# ======================
# Build /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/chr_marker/grape_chrMarkers_allele_db_v001.fa
# ======================


# =====================
# Concatenate chr and nonChr markers into one file
# =====================
```bash
cat /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/nonChr_marker_coreDB/grape_nonChr_allele_db_v001.fa \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/chr_marker_coreDB/grape_chrMarkers_allele_db_v001.fa \
> /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/grape_allele_db_v001.fa
```


# Concatenate matchCnt_lut.txt



# Test code
```bash
py /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db05_determine_alleleOri_from_blast_AND_update_f180bp_v1.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10633_Vvinifera/OrderAppendix_1_DG25-10633_flipRefAlt_snpID_ref_alt_amplicons.fa.f300bp.bn \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys_ref_alt.fa
```


# Prepare vcf header file




# Find all unique paralogous alleles between files
```bash
python3 /Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/04_util/util_concat_2files_keepUnique.py \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/db_versions/grape_allele_db_v002_DG25-10632_paralogousHaps.csv \
/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/db_versions/grape_allele_db_v003_DG25-10633_paralogousHaps.csv
```
[INFO] Reading first allele file: /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/db_versions/grape_allele_db_v002_DG25-10632_paralogousHaps.csv
[INFO] Reading second allele file: /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/db_versions/grape_allele_db_v003_DG25-10633_paralogousHaps.csv
[INFO] Combining rows (File 1: 967 rows | File 2: 1,301 rows)
[INFO] Parsing AlleleID prefixes before the ':' delimiter for unique targeting...
[INFO] Deduplicating dataset based on extracted sub-string unique identifier...
[INFO] Master consolidation finalized: 1,661 unique alleles preserved.
[INFO] Exporting final CSV matrix to: /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/db_versions/grape_allele_db_v002_DG25-10632_paralogousHaps_combined_alleles.csv
[INFO] Processing complete.


#


