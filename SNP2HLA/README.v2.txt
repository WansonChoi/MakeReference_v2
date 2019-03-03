#################################################################################################

    < SNP2HLA.py >

    SNP2HLA: Imputation of HLA amino acids and classical alleles from SNP genotypes

    Author: Sherman Jia (xiaomingjia@gmail.com)
            + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
            + Extensive modifications by Phil Stuart (pstuart@umich.edu) on 4/19/16 to allow use of Beagle 4.1:
                verified to work with 22Feb16, 15Apr16, and 03May16 versions of Beagle.
            + Small modifications by Yang Luo (yangluo@broadinstitute.org): 09/30/16: verfiied working with Bealge 4.1 27Jun16.
            + Recoded to Python and updated by Wanson Choi(wschoi.bhlab@gmail.com) : 2019/02/06 
            
            
    DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.        
    
    INPUTS:
    1. Plink dataset (*.bed/bim/fam)
    2. Reference dataset (*.bgl.phased.vcf.gz(Beagle 4.1), *.markers(Beagle 3.0.4); *.fam/.bim/.FRQ.frq in PLINK format)

    DEPENDENCIES: (download and place in the same folder as this script)
    1. PLINK (1.9)  (Will not work with older Plink 1.07)
    2. Beagle (4.1) (Need to rename java executable as beagle.jar)
    3. vcf2gprobs.jar (Beagle utility for generating a Beagle v3 genotypes probability file from a Beagle 4.1 vcf file with GT field data)
    4. [Optional] If genetic_map_file argument is specified, PLINK format genetic map on cM scale (plink.chr6.GRCh36.map, downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)

    cf) 
    - 'revert_alleles' is no longer needed because the problems related to unacceptable allele characters and duplicated genomic positions are preprocessed in MakeReference_v2.
    - 'beagle2vcf.jar' is also no longer needed in SNP2HLA because the conversion from *.bgl to *.vcf is done in MakeReference_v2.
    - 'merge_tables.pl' and 'ParseDosage.csh' sources are prepared in 'src/' directory of project folder.

    
    USAGE:
    python3 SNP2HLA.py \
            --input `DATA (.bed/.bim/.fam)` \
            --reference `REFERENCE (.bgl.phased.vcf.gz/.markers/.fam/.bim/.bed/.FRQ.frq)` \
            --out `OUTPUT`

    (ex1)
    python3 SNP2HLA.py \
            --input data/1958BC \
            --reference data/Reference_Panel_bglv4/HM_CEU_REF.hg18.imgt3320.bglv4 \
            --out TEST_ex1_SNP2HLA

    (ex2)
    python3 SNP2HLA.py \
            --input data/1958BC \
            --reference data/Reference_Panel_bglv4/HM_CEU_REF.hg18.imgt3320.bglv4 \
            --out TEST_ex2_SNP2HLA \
            --dependency dependency/ \
            --java-mem 10G \
            --nthreads 4 \
            --iter 10

    To check more detailed arguments, Type "python3 SNP2HLA.py --help"


#################################################################################################



Thank you for downloading SNP2HLA. To use this package:

1. Prepare dependent software.

    (1-1) Plink 1.9
    - https://www.cog-genomics.org/plink2/

    (1-2) Beagle 4.1 
    - https://faculty.washington.edu/browning/beagle/b4_1.html

    (1-3) Beagle Utilities (vcf2gprobs.jar)
    - https://faculty.washington.edu/browning/beagle_utilities/utilities.html#vcf2gprobs


    - Download those 3 software which are proper to your system(ex. Linux or OS X).
    - In case of Beagle4.1, You must rename file name of downloaded beagle 4.1 jar file to "beagle.jar". 
      (i.e. "beagle.27Jan18.7e1.jar" -> "beagle.jar").
    - Copy and Paste them into either (1) SNP2HLA project folder or (2) Somewhere directory that you'd like to manage those software.
    - If you choose (1) way, then don't pass "--dependency" option, It will automatically recognize project folder as where dependent software are.
    - If you choose (2) way, then pass the path of that directory to "--dependency", then it will search dependent software in there.



2. Run SNP2HLA with sample data provided (10 samples from British 1958 Birth Cohort, build 35, HapMap CEU reference dataset) using the following command:

    python3 SNP2HLA.py \
            --input data/1958BC \
            --reference data/Reference_Panel_bglv4/HM_CEU_REF.hg18.imgt3320.bglv4 \
            --out TEST_ex2_SNP2HLA \
            --dependency dependency/ \
            --java-mem 10G \
            --marker-window 1000 \
            --nthreads 4 \
            --iter 10

In the above example,
- ("--input") 1958BC is the SNP genotype plink files (.bed/.bim/.fam),
- ("--reference") HM_CEU_REF.hg18.imgt3320.bglv4 is the reference dataset (.bgl.phased.vcf.gz/.markers)
- ("--dependency") "dependency/" folder would be the directory where the user gathered dependent software.
- ("--java-mem") 10G(gigabyte) is the maximum java heap size (in mb) for imputation using Beagle (increase as needed)
- ("--marker-window") 1000 is the marker window size that Beagle uses for phasing and imputation
- ("--nthreads") The number of threads that will be used in Beagle4.1 imputation.
- ("--iter") The number of iteration which will be done in imputation.


SNP2HLA will also run with default parameters if memory and window size are not provided:
(java mamory = 2Gb, marker window size = 1000)

----------------------------------------------------------------------------------

Files included in this package:

1. SNP2HLA.py: Performs imputation (via Beagle4.1) after SNP QC (using PLINK)
2. Merge_tables.pl: Merges files according to indices in a particular column (called by SNP2HLA.py)
3. ParseDosage.csh: Converts .gprobs (Beagle) file to .dos (PLINK) file
4. HapMap CEU reference dataset (Plink and Beagle formats)
5. Sample SNP dataset of 10 individuals from Britist 1958 Birth Cohort (1958BC.bed/.bim/.fam)

----------------------------------------------------------------------------------

Input files (for SNP2HLA.py):

1. SNP dataset (.bed/bim/fam PLINK format)
   *** We compare rsIDs to Reference, so coordinates (hg18/hg19) are not important. ***
2. Reference dataset (.bgl.phased.vcf.gz/.markers Beagle format)
3. Pointer to dependent software(plink1.9, beagle4.1, vcf2gprobs.jar)


Output files:

- {OUTPUT}.dosage: PLINK format dosage data (recommended for downstream analysis)
- {OUTPUT}[.bed/.bim/.fam]: PLINK format best-guess genotype files
- {OUTPUT}.bgl.gprobs: imputation posterior probabilities for SNPs, HLA alleles, and HLA amino acids
- {OUTPUT}.bgl.phased.vcf.gz: imputation best-guess genotypes 
*** Output coordinates are all in hg18 currently. ***

cf) Information of {OUTPUT}.bgl.r2 file which was created in original version of SNP2HLA is now 
    included in {OUTPUT}.bgl.phased.vcf.gz as "AR2" of "INFO" column.

----------------------------------------------------------------------------------

Pivotal rule for using allele characters in MakeReference_v2 and SNP2HLA is same as the one of original MakeReference and SNP2HLA.

    - Marker Nomenclature: For binary encodings, P = Present, A = Absent.

However, as Beagle v3.0.1 is replaced by Beagle v4.1, Two major issues come out and must be solved.

    (1) No other allele characters execpt 'A', 'C', 'G', 'T' and 'N' are allowed,
    (2) No same genomic position is allowed to each markers.

To solve these issues in using Beagle v4.1, MakeReference_v2 does some temporary preprocessing to its output reference panel.

    (1) Any bi-allelic characters which are not composed of 'A', 'C', 'G', 'T' and 'N' will be encoded to {'A', 'T'} temporarily.
    (2) Markers which have same genomic position will get a new temporary genomic position which is unique and nearest to its original position.

For example, if {'P', 'A'} is given, then it will be encoded to {'T', 'A'}. Likewise, if {'L', 'V'} is given, then it will be encoded to {'T', 'A'}. 
Those encoding information which allele characters are encoded to {'A', 'T'} will be stored in "*.pENCODED.aENCODED.emap", one file of reference panel. 
Users will be able to revert them back or interpret them based on it.

In summary, even though it looks there are no more {'P', 'A'} bi-allelic characters, it could have been encoded to {'T', 'A'} and the old rule 
using 'P' as 'Present' and 'A' as 'Absent' is still valid.




1. Classical HLA alleles: HLA_[GENE]*[ALLELE]. 
- HLA_C*03:04:01:01 = HLA-C*03:04:01:01 (Standard 4-field allele; Here, 'standard' means it is declared by IMGT-HLA organization.)
- HLA_DRB1*07 = HLA-DRB1*07 (1-field allele (2 to 3 digit allele))

2. HLA Amino Acids: AA_[GENE]_[AMINO ACID RELATIVE POSITION]_[GENETIC POSITION]_[Nth EXON or INTRON]_[ALLELE]. 
- AA_A_9_30018537_exon2_F = 9th amino acid of 2nd exon of HLA-A, genetic position 30018537 (center of codon), allele = F (Phe) of multi-allelic position
- AA_C_291_31345800_exon5 = 291th amino acid 5th exon of HLA-C, genetic position 31345800, bi-allelic 
    (check both "HM_CEU_REF.hg18.imgt3320.pENCODED.aENCODED.emap" and "{OUTPUT}.bim" for alleles, P = Present, A = Absent)

3. HLA intragenic SNPS: SNPS_[GENE]_[DNA SEQUENCE RELATIVE POSITION]_[GENOMIC POSITION]_[Nth EXON or INTRON]_[ALLELE]
- SNPS_B_2596_31430319_intron6_G = 2596th SNP at position 31430319 of 6th intron of HLA-B, allele = G (guanine) of multi-allelic position
- SNPS_DRB1_5499_32659999_exon2 = = 5499th SNP at position 32659999 of 2nd exon of HLA-DRB1, bi-allelic 
    (check both "HM_CEU_REF.hg18.imgt3320.pENCODED.aENCODED.emap" and "{OUTPUT}.bim" for alleles, P = Present, A = Absent)
- SNPS_DQB1_1436_32740927_intron1_AC = 1436th SNP at position 32740927 of 1st intron of HLA-DQB1, alleles = A (adenine) or C (Cysteine), 
    (check both "HM_CEU_REF.hg18.imgt3320.pENCODED.aENCODED.emap" and "{OUTPUT}.bim" for alleles, P = Present, A = Absent)

4. INDEL(Insertions / deletions): INDEL_[VARIANT(AA/SNPS)]_[GENE]_[RELATIVE POSITION]_[GENOMIC POSITION]
- Z : Insertion / z : Deletion in "HM_CEU_REF.hg18.imgt3320.pENCODED.aENCODED.emap" (anyway they will appear encoded to {'T', 'A'} in "{OUTPUT}.bim")
- INDEL_AA_C_300x301_31345772 = Indel at amino acids between relative position 300 and 301 in HLA-C, genetic position 31345772 
    (check both "HM_CEU_REF.hg18.imgt3320.pENCODED.aENCODED.emap" and "{OUTPUT}.bim" for alleles, Z = Insertion, z = Deletion)
- INDEL_SNPS_C_2033x2034_31345794 = Indel at intragenic SNPs between relative position 2033 and 2034 of HLA-C, genetic position 31345794
    (check both "HM_CEU_REF.hg18.imgt3320.pENCODED.aENCODED.emap" and "{OUTPUT}.bim" for alleles, Z = Insertion, z = Deletion)
- INDEL_SNPS_DQB1_5365x5366_32736998 = Indel at amino acids between relative position 5365 and 5366 of HLA-DQB1, genetic position 32736998
    (check both "HM_CEU_REF.hg18.imgt3320.pENCODED.aENCODED.emap" and "{OUTPUT}.bim" for alleles, Z = Insertion, z = Deletion)

----------------------------------------------------------------------------------

Association testing in PLINK

plink --noweb --dosage OUTPUT.dosage noheader format=1 --fam OUTPUT.fam --logistic --out OUTPUT.assoc

----------------------------------------------------------------------------------
