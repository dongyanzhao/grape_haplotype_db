#!/bin/bash
# For 81-bp amplicons from DArTag reports
#  1. Created the 180 bp flanking sequence db and update sequences orientation based on DArTag Ref and Alt amplicons BLAST result
   #    1). Prepare 180bp flanking sequences of Ref and Alt alleles from probe design file submitted to DArT
   #    2). Extract amplicon sequences of Ref and Alt alleles from DArTag report and generate MATCH allele LUT
   #    3). BLAST Ref and Alt amplicon sequences from DArTag report to the 180 bp flanking sequences
   #    4). Determine alignment orientation between amplicon vs. f180-bp, based on which update the f180-bp sequences
   #    5). Rerun BLAST against the updated flankseq sequences
   #    6). Get sfetch keys for the amplicons (because the 3' end of some amplicons are inaccurate)


# change here ######
exec &> /Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/chr_marker_coreDB/grape_chrMarkers_allele_db_v001_process.readme

# change here ######
SCRIPTS_DIR='/Users/dongyan.zhao/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB'
MARKERID_LUT='/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut.csv'
ALLELE_DB_DIR='/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/chr_marker_coreDB'
ALLELE_DB='grape_chrMarkers_allele_db_v001.fa'
MATCH_CNT='grape_allele_db_v001_matchCnt_lut.txt'
REPORT='/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_P00_validation/data/Report-DG25-10632_Margaret_6plates/DG25-10632_MADC_flipRefAlt.csv'
FLANK_SEQ='/Users/dongyan.zhao/PycharmProjects/BI/grape_dartag_00_microhaplotype_db/data/flankseq/20250904-Grape_3K_DArTag_probe_rm4Dup_hyphenAnno_rm3nonMatch_rm5miss_chrMarkers_lut_f300bp_sfetchKeys_ref_alt.fa'
FLANK_LEN='f300bp'
REF_LEN=250

now=$(date)
printf "%s\n" "$now"
printf '\n# 1). Prepare 180bp flanking sequences of Ref and Alt alleles from Marker IDs'
printf '\n  # Processed separately because there are IUPAC codes in the probe design file'
printf '\n  # If there is already a SNP ID LUT, provide "N" in the command line\n'
#python3 $SCRIPTS_DIR/db01_get_f180bp_fasta_AND_snpID_lut_from_probeDesign_v1.py $PROBE N
#FLANK_SEQ=${PROBE%????}'_f180bp.fa'


printf  '\n# 2). Update snpIDs in DArTag report\n'
python3 $SCRIPTS_DIR/db02_update_snpID_in_madc_v1.py $MARKERID_LUT $REPORT
REPORT_ID=${REPORT%????}'_snpID.csv'


printf  '\n# 3). Extract amplicon sequences of Ref and Alt alleles from DArTag report and generate MATCH allele LUT\n'
python3 $SCRIPTS_DIR/db03_ext_ref_alt_amp_AND_gen_match_lut_from_madc_v1.py $REPORT_ID
REF_ALT=${REPORT%????}'_snpID_ref_alt_amplicons.fa'


printf  '\n# 4). BLAST Ref and Alt amplicon sequences from DArTag report to the xx bp flanking sequences\n'
#makeblastdb -in $FLANK_SEQ -dbtype nucl
# 20250404_grape_3K_snpID_lut_f300bp_sfetchKeys_ref_alt.fa
REF_ALT_BLAST=${REPORT%????}'_snpID_ref_alt_amplicons.fa.'$FLANK_LEN'.bn'
MATCHCNT_LUT=${REPORT%????}'_snpID_matchCnt_lut.txt'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ -query $REF_ALT -out $REF_ALT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 10 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
### 2025.3.21 I had to change the max_target_seqs to 10 because some marker are close to each other and the amplicons are similar 


printf  '\n# 5). Determine alignment orientation between amplicon vs. f180-bp, based on which update the f180-bp sequences\n'
python3 $SCRIPTS_DIR/db05_determine_alleleOri_from_blast_AND_update_f180bp_v1.py $REF_ALT_BLAST $FLANK_SEQ


printf  '\n# 6). Rerun BLAST against the updated flankseq sequences\n'
FLANK_SEQ_REV=${FLANK_SEQ%???}'_rev.fa'
makeblastdb -in $FLANK_SEQ_REV -dbtype nucl
REF_ALT_BLAST_REV=${REPORT%????}'_snpID_ref_alt_amplicons.fa.'$FLANK_LEN'_rev.bn'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ_REV -query $REF_ALT -out $REF_ALT_BLAST_REV -evalue 1e-5 -num_threads 6 -max_target_seqs 10 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


printf  '\n# 7). Get sfetch keys for the amplicons (because the 3 end of some amplicons are inaccurate)\n'
python3 $SCRIPTS_DIR/db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py $MARKERID_LUT $REF_ALT_BLAST_REV $REF_LEN


printf  '\n# 8). Get the amplicon sequences from the flanking sequences\n'
esl-sfetch --index $FLANK_SEQ_REV
REF_ALT_BLAST_REV_SFETCH=${REPORT%????}'_snpID_ref_alt_amplicons.fa.'$FLANK_LEN'_rev.bn_'$REF_LEN'bp_sfetchKeys.txt'
SFETCH_FA=${REPORT%????}'_snpID_ref_alt_amplicons.fa.'$FLANK_LEN'_rev_'$REF_LEN'bp_sfetch.fa'
esl-sfetch -Cf $FLANK_SEQ_REV $REF_ALT_BLAST_REV_SFETCH > $SFETCH_FA


printf  '\n# 9) Copy allele fasta and match count lut to db directory\n'
cp $SFETCH_FA $ALLELE_DB_DIR/$ALLELE_DB
cp $MATCHCNT_LUT $ALLELE_DB_DIR/$MATCH_CNT


printf  '\n# 10) Make blastdb\n'
makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB -dbtype nucl
