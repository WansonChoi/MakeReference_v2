# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
from pathlib import Path
from platform import platform

########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def MakeReference(_INPUT_DATA, _hped, _OUTPUT,
                  _dictionary_AA, _dictionary_SNPS,
                  _previous_version=False, _hg="19",
                  _mem="2000m", _p_depedency="./dependency",
                  __save_intermediates=False):

    """
    """

    """
    (2019. 1. 5.) by Wanson
    < Main arguments of MakeReference modification >
    - `_p_plink`, `_p_beagle`, `_p_linkage2beagle` (dependency related arguments) 
        => integrated to `_p_dependency`.
    - `_dictionary_AA_map` and `_dictionary_SNPS_map` 
        => integrated to `_dictionary_AA` and `_dictionary_SNPS`. Each two dictionary files will be dealt with the file prefix
        (ex. "data/HLA_DICTIONARY_AA.txt" and "data/HLA_DICTIONARY_AA.map" => dealt with "data/HLA_DICTIONARY_AA"
        => now it is positional arguments not optional.
        (ex. "data/HLA_DICTOINARY_SNPS.txt" and "data/HLA_DICTOINARY_SNPS.map" => dealt with "data/HLA_DICTIONARY_SNPS"
    - `_previous_version` set "True" by default.
    - `_hg` set to "18" by default.    
    
    """

    ########## < Preprocessing main arguments > ##########

    if _previous_version:

        # Human Genome version.
        if _hg != "18":
            _hg = "18"

        # Dictionary (Use default dictionary)
        if not (_dictionary_AA and _dictionary_SNPS):
            _dictionary_AA = "data/MakeReference_old/HLA_DICTIONARY_AA"
            _dictionary_SNPS = "data/MakeReference_old/HLA_DICTIONARY_SNPS"



    # Memory representation
    p = re.compile(r'\d+m$')
    p2 = re.compile(r'g$')

    if not p.match(_mem):

        if p2.search(_mem):
            _mem = p.sub(repl="000m", string=_mem) # Gigabyte to Megabyte to use it in java.
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Given memory for beagle is unappropriate. Please check it again.\n")
            sys.exit()


    ########## < Path Variables > ##########

    ### Major path variables for project

    # src
    p_src_MakeReferece = "./src/MakeReference_old" if _previous_version else "./src/MakeReference"


    ########## < CHECK FOR DEPENDENCIES > ##########

    ### Other Software.

    _p_plink = os.path.join(_p_depedency, "plink_mac" if not bool(re.search(pattern="Linux", string=platform())) else "plink") #plink v1.9
    _p_beagle = os.path.join(_p_depedency, "beagle.jar")
    _p_linkage2beagle = os.path.join(_p_depedency, "linkage2beagle.jar")
    _p_beagle2vcf = os.path.join(_p_depedency, "beagle2vcf")

    if not os.path.exists(_p_plink):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'PLINK' (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) in '{}'\n".format(_p_depedency))
        sys.exit()
    if not os.path.exists(_p_beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'Beagle 4' (https://faculty.washington.edu/browning/beagle/b4_0.html#download) in '{}'\n".format(_p_depedency)) # Now we use beagle 4.
        sys.exit()
    if not os.path.exists(_p_linkage2beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'linkage2beagle.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) in '{0}'\n".format(_p_depedency))
        sys.exit()
    if not os.path.exists(_p_beagle2vcf):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please Prepare 'beagle2vcf.jar' (http://faculty.washington.edu/browning/beagle_utilities/utilities.html) in '{0}'\n".format(_p_depedency))
        sys.exit()


    ### Dictionary Information for HLA sequence

    _dictionary_AA_seq = _dictionary_AA + ".txt"
    _dictionary_AA_map = _dictionary_AA + ".map"

    _dictionary_SNPS_seq = _dictionary_SNPS + ".txt"
    _dictionary_SNPS_map = _dictionary_SNPS + ".map"

    if not os.path.exists(_dictionary_AA_map):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please check again whether '{}' is given properly.\n".format(_dictionary_AA_map))
        sys.exit()
    if not os.path.exists(_dictionary_AA_seq):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please check again whether '{}' is given properly.\n".format(_dictionary_AA_seq))
        sys.exit()
    if not os.path.exists(_dictionary_SNPS_map):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please check again whether '{}' is given properly.\n".format(_dictionary_SNPS_map))
        sys.exit()
    if not os.path.exists(_dictionary_SNPS_seq):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please check again whether '{}' is given properly.\n".format(_dictionary_SNPS_seq))
        sys.exit()


    ### Source Code Scripts

    if not _previous_version:

        # New version with Python.

        if not os.path.exists(os.path.join(p_src_MakeReferece, "HLAtoSequences.py")):
            print(std_ERROR_MAIN_PROCESS_NAME + "Error. 'HLAtoSequences.py' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()
        else:
            from src.MakeReference.HLAtoSequences import HLAtoSequences

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeVariants.py")):
            print(std_ERROR_MAIN_PROCESS_NAME + "Error. 'encodeVariants.py' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()
        else:
            from src.MakeReference.encodeVariants import encodeVariants

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeHLA.py")):
            print(std_ERROR_MAIN_PROCESS_NAME + "Error. 'encodeHLA.py' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()
        else:
            from src.MakeReference.encodeHLA import encodeHLA

    else:
        # Previous version with Perl.

        if not os.path.exists(os.path.join(p_src_MakeReferece, "HLAtoSequences.pl")):
            print(std_ERROR_MAIN_PROCESS_NAME + "Error. 'HLAtoSequences.pl' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeVariants.pl")):
            print(std_ERROR_MAIN_PROCESS_NAME + "Error. 'encodeVariants.pl' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()

        if not os.path.exists(os.path.join(p_src_MakeReferece, "encodeHLA.pl")):
            print(std_ERROR_MAIN_PROCESS_NAME + "Error. 'encodeHLA.pl' not found in '{0}'".format(p_src_MakeReferece))
            sys.exit()



    ########## <Core Variables> ##########

    HLA_DATA = _hped


    ### Intermediate path.
    OUTPUT = _OUTPUT if not _OUTPUT.endswith('/') else _OUTPUT.rstrip('/')
    if bool(os.path.dirname(OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        # If `os.path.dirname(OUTPUT)` doesn't exist, then it means the output of MakeReference should be genrated in current directory.
        INTERMEDIATE_PATH = "./"


    # (2018. 5. 29) additional processing for intermediate path and output file prefix
    SNP_DATA = _INPUT_DATA
    SNP_DATA2 = os.path.join(INTERMEDIATE_PATH, os.path.basename(SNP_DATA)) # for each output directory


    plink = ' '.join([_p_plink, "--allow-no-sex","--silent"])
    beagle = ' '.join(["java", " -Xmx"+_mem, "-Xss5M", "-jar", _p_beagle])
    linkage2beagle = ' '.join(["java", "-Xmx"+_mem, "-jar", _p_linkage2beagle])
    beagle2vcf = ' '.join(["java", "-Xmx"+_mem, "-jar", _p_beagle2vcf])


    ########## <Flags for Code Block> ##########

    ENCODE_AA = 1
    ENCODE_HLA = 0
    ENCODE_SNPS = 0

    EXTRACT_FOUNDERS = 0
    MERGE = 0
    QC = 0

    PREPARE = 0
    PHASE = 0
    CLEANUP = 0 # set to zero for time being



    ########## <Making Reference Panel> ##########

    print(std_MAIN_PROCESS_NAME + "Making Reference Panel for \"{0}\"".format(OUTPUT))

    index = 1

    if ENCODE_AA:

        '''
        echo "[$i] Generating amino acid sequences from HLA types.";  @ i++
        ./HLAtoSequences.pl $HLA_DATA HLA_DICTIONARY_AA.txt AA > $OUTPUT.AA.ped
        cp HLA_DICTIONARY_AA_hg19.map $OUTPUT.AA.map # hg19
        # cp HLA_DICTIONARY_AA.map $OUTPUT.AA.map

        echo "[$i] Encoding amino acids positions." ;  @ i++
        ./encodeVariants.pl $OUTPUT.AA.ped $OUTPUT.AA.map $OUTPUT.AA.CODED

        plink --file $OUTPUT.AA.CODED --missing-genotype 0 --make-bed --out $OUTPUT.AA.TMP
        awk '{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}' $OUTPUT.AA.TMP.bim | grep -v INS | cut -f2 > to_remove
        plink --bfile $OUTPUT.AA.TMP --exclude to_remove --make-bed --out $OUTPUT.AA.CODED

        # rm $OUTPUT.AA.TMP*; rm to_remove
        # rm $OUTPUT.AA.???
        '''

        print("[{}] Generating amino acid sequences from HLA types.".format(index))

        if not _previous_version:
            HLAtoSequences(HLA_DATA, _dictionary_AA_seq, "AA", _out=OUTPUT)
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "HLAtoSequences.pl"), HLA_DATA, _dictionary_AA, "AA", ">", OUTPUT+".AA.ped"])
            print(command)
            os.system(command)

        os.system(' '.join(["cp", _dictionary_AA_map, OUTPUT + '.AA.map']))

        index += 1



        print("[{}] Encoding amino acids positions.".format(index))

        if not _previous_version:
            encodeVariants(OUTPUT + '.AA.ped', OUTPUT + '.AA.map', OUTPUT + '.AA.CODED') # previously "enCODED".
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "encodeVariants.pl"), OUTPUT+".AA.ped", OUTPUT+".AA.map", OUTPUT+".AA.CODED"])
            print(command)
            os.system(command)

        # command for checking output from encodeVariant.py(.pl)
        command = ' '.join([plink, "--file", OUTPUT+'.AA.CODED', "--missing-genotype 0", "--make-bed", "--out", OUTPUT+'.AA.TMP'])
        print(command)
        os.system(command)

        command = ' '.join(["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT + '.AA.TMP.bim', "|", "grep -v {0}".format("INDEL" if not _previous_version else "INS"), "|", "cut -f2", ">", os.path.join(INTERMEDIATE_PATH, "to_remove")])

        """
        In previous framework originally created by Sherman Jia, Only insertions were dealt with as a marker "INS".
        In the new version of Framework, marker label is "INDEL".
        """

        print(command)
        os.system(command)

        command = ' '.join([plink, "--bfile", OUTPUT+'.AA.TMP', "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"), "--make-bed", "--out", OUTPUT+'.AA.CODED'])
        print(command)
        os.system(command)

        index += 1


        if not __save_intermediates:

            os.system("rm " + (OUTPUT + ".AA.{ped,map}"))
            os.system("rm " + (OUTPUT + ".AA.TMP.*"))
            os.system("rm " + os.path.join(INTERMEDIATE_PATH, "to_remove"))
            os.system("rm " + OUTPUT + ".AA.CODED.{ped,map}")




    if ENCODE_HLA:

        print("[{}] Encoding HLA alleles.".format(index))

        if not _previous_version:
            encodeHLA(HLA_DATA, OUTPUT, _hg)
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "encodeHLA.pl"), HLA_DATA, OUTPUT+".HLA.map", ">", OUTPUT+".HLA.ped"])
            print(command)
            os.system(command)


        command = ' '.join([plink, "--file", OUTPUT+'.HLA', "--make-bed", "--out", OUTPUT+'.HLA'])
        print(command)
        os.system(command)

        index += 1



    if ENCODE_SNPS:

        print("[{}] Generating DNA sequences from HLA types.".format(index))

        if not _previous_version:
            HLAtoSequences(HLA_DATA, _dictionary_SNPS_seq, "SNPS", OUTPUT)
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "HLAtoSequences.pl"), HLA_DATA, _dictionary_SNPS, "SNPS", ">", OUTPUT+".SNPS.ped"])
            print(command)
            os.system(command)

        command = ' '.join(["cp", _dictionary_SNPS_map, OUTPUT+'.SNPS.map'])
        print(command)
        os.system(command)

        index += 1


        print("[{}] Encoding SNP positions.".format(index))

        if not _previous_version:
            encodeVariants(OUTPUT+'.SNPS.ped', OUTPUT+'.SNPS.map', OUTPUT+'.SNPS.CODED')
        else:
            command = ' '.join([os.path.join(p_src_MakeReferece, "encodeVariants.pl"), OUTPUT+".SNPS.ped", OUTPUT+".SNPS.map", OUTPUT+".SNPS.CODED"])
            print(command)
            os.system(command)


        command = ' '.join([plink, "--file", OUTPUT+'.SNPS.CODED', "--missing-genotype 0", "--make-bed", "--out", OUTPUT+'.SNPS.TMP'])
        print(command)
        os.system(command)

        command = ' '.join(["awk", '\'{if ($5 == "0" || $5 == "x" || $6 == "x"){print $2}}\'', OUTPUT +'.SNPS.TMP.bim', "|", "grep -v {0}".format("INDEL" if not _previous_version else "INS"), "|", "cut -f2", ">", os.path.join(INTERMEDIATE_PATH, "to_remove")])
        print(command)
        os.system(command)

        command = ' '.join([plink, "--bfile", OUTPUT+'.SNPS.TMP', "--exclude", os.path.join(INTERMEDIATE_PATH, "to_remove"), "--make-bed", "--out", OUTPUT+'.SNPS.CODED'])
        print(command)
        os.system(command)

        index += 1


        rm_tlist = (OUTPUT+'.SNPS.TMP*', os.path.join(INTERMEDIATE_PATH, 'to_remove'), OUTPUT+'.SNPS.???')

        for i in rm_tlist:
            print(i)
            os.system("rm "+i)



    if EXTRACT_FOUNDERS:

        print("[{}] Encoding SNP positions.".format(index))

        """
        if ($EXTRACT_FOUNDERS) then
            echo "[$i] Extracting founders."; @ i++
            # founder의 정의가 이게 맞는건 아니겠지만, 아래 plink명령어를 거치고 나오는 founder라는 애들은 모두 엄마,아빠 ID정보가 없는 애들임.
            plink --bfile $SNP_DATA --filter-founders --mind 0.3 --alleleACGT --make-bed --out $SNP_DATA.FOUNDERS

            # Initial QC on Reference SNP panel
            plink --bfile $SNP_DATA.FOUNDERS --hardy        --out $SNP_DATA.FOUNDERS.hardy  # 진짜 92명에 대해 position별로 HWE test한 결과
            plink --bfile $SNP_DATA.FOUNDERS --freq         --out $SNP_DATA.FOUNDERS.freq   # 실제 --freq 옵션이 allele frequency계산해주는 옵션임.
            plink --bfile $SNP_DATA.FOUNDERS --missing      --out $SNP_DATA.FOUNDERS.missing
            awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.hardy.hwe      | awk ' $9 < 0.000001 { print $2 }' | sort -u > remove.snps.hardy
            awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.freq.frq       | awk ' $5 < 0.01 { print $2 } '             > remove.snps.freq
            awk '{if (NR > 1){print}}' $SNP_DATA.FOUNDERS.missing.lmiss  | awk ' $5 > 0.05 { print $2 } '              > remove.snps.missing
            cat remove.snps.*                                            | sort -u                                     > all.remove.snps

            plink --bfile $SNP_DATA.FOUNDERS --allow-no-sex --exclude all.remove.snps --make-bed --out $SNP_DATA.FOUNDERS.QC

            # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

            plink --bfile $OUTPUT.HLA --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.HLA.FOUNDERS
            plink --bfile $OUTPUT.SNPS.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.SNPS.FOUNDERS
            plink --bfile $OUTPUT.AA.CODED --filter-founders --maf 0.0001 --make-bed --out $OUTPUT.AA.FOUNDERS

            rm remove.snps.*
        endif
        """
        SNP_DATA = _INPUT_DATA
        _INPUT_DATA_prefix = Path(_INPUT_DATA).name
        command = ' '.join([plink, "--bfile", SNP_DATA, "--filter-founders", "--mind 0.3", "--alleleACGT", "--make-bed", "--out", os.path.join(INTERMEDIATE_PATH, _INPUT_DATA_prefix+'.FOUNDERS')])
        print(command)
        os.system(command)

        # SNP_DATA2 = os.path.join(INTERMEDIATE_PATH, _INPUT_DATA_prefix)

        # Initial QC on Reference SNP panel
        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--hardy", "--out", SNP_DATA2+'.FOUNDERS.hardy'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--freq", "--out", SNP_DATA2+'.FOUNDERS.freq'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--missing", "--out", SNP_DATA2+'.FOUNDERS.missing'])
        print(command)
        os.system(command)

        command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.hardy.hwe', "|", "awk", "' $9 < 0.000001 { print $2 }'", "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.hardy")])
        print(command)
        os.system(command)
        command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.freq.frq', "|", "awk", "' $5 < 0.01 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.freq")])
        print(command)
        os.system(command)
        command = ' '.join(["awk", "'{if (NR > 1){print}}'", SNP_DATA2+'.FOUNDERS.missing.lmiss', "|", "awk", "' $5 > 0.05 { print $2 } '", ">", os.path.join(INTERMEDIATE_PATH, "remove.snps.missing")])
        print(command)
        os.system(command)
        command = ' '.join(["cat", os.path.join(INTERMEDIATE_PATH, "remove.snps.*"), "|", "sort -u", ">", os.path.join(INTERMEDIATE_PATH, "all.remove.snps")])
        print(command)
        os.system(command)


        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS', "--exclude", os.path.join(INTERMEDIATE_PATH, "all.remove.snps"), "--make-bed", "--out", SNP_DATA2+'.FOUNDERS.QC'])
        print(command)
        os.system(command)

        # Founders are identified here as individuals with "0"s in mother and father IDs in .fam file

        command = ' '.join([plink, "--bfile", OUTPUT+'.HLA', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.HLA.FOUNDERS'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", OUTPUT+'.SNPS.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.SNPS.FOUNDERS'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", OUTPUT+'.AA.CODED', "--filter-founders", "--maf 0.0001", "--make-bed", "--out", OUTPUT+'.AA.FOUNDERS'])
        print(command)
        os.system(command)

        index += 1


        command = ' '.join(["rm", os.path.join(INTERMEDIATE_PATH, "remove.snps.*")])
        print(command)
        os.system(command)



    if MERGE:

        print("[{}] Merging SNP, HLA, and amino acid datasets.".format(index))

        """
        echo "[$i] Merging SNP, HLA, and amino acid datasets.";  @ i++
        echo "$OUTPUT.HLA.FOUNDERS.bed $OUTPUT.HLA.FOUNDERS.bim $OUTPUT.HLA.FOUNDERS.fam" > merge_list
        echo "$OUTPUT.AA.FOUNDERS.bed $OUTPUT.AA.FOUNDERS.bim $OUTPUT.AA.FOUNDERS.fam" >> merge_list
        echo "$OUTPUT.SNPS.FOUNDERS.bed $OUTPUT.SNPS.FOUNDERS.bim $OUTPUT.SNPS.FOUNDERS.fam" >> merge_list
        plink --bfile $SNP_DATA.FOUNDERS.QC --merge-list merge_list --make-bed --out $OUTPUT.MERGED.FOUNDERS
        rm $OUTPUT.HLA.???
        rm $OUTPUT.AA.CODED.???
        rm $OUTPUT.SNPS.CODED.???
        rm merge_list

        """

        TMP_merged_list = os.path.join(INTERMEDIATE_PATH, "merge_list")

        command = ' '.join(["echo", OUTPUT+'.HLA.FOUNDERS.bed', OUTPUT+'.HLA.FOUNDERS.bim', OUTPUT+'.HLA.FOUNDERS.fam', ">", TMP_merged_list])
        print(command)
        os.system(command)

        command = ' '.join(["echo", OUTPUT+'.AA.FOUNDERS.bed', OUTPUT+'.AA.FOUNDERS.bim', OUTPUT+'.AA.FOUNDERS.fam', ">>", TMP_merged_list])
        print(command)
        os.system(command)

        command = ' '.join(["echo", OUTPUT+'.SNPS.FOUNDERS.bed', OUTPUT+'.SNPS.FOUNDERS.bim', OUTPUT+'.SNPS.FOUNDERS.fam', ">>", TMP_merged_list])
        print(command)
        os.system(command)

        command = ' '.join([plink, "--bfile", SNP_DATA2+'.FOUNDERS.QC', "--merge-list", TMP_merged_list, "--make-bed", "--out", OUTPUT+'.MERGED.FOUNDERS'])
        print(command)
        os.system(command)

        index += 1


        rm_tlist = (OUTPUT+'.HLA.???', OUTPUT+'.AA.CODED.???', OUTPUT+'.SNPS.CODED.???', TMP_merged_list)

        for i in rm_tlist:
            print(i)
            os.system("rm "+i)



    if QC:

        print("[{}] Performing quality control.".format(index))

        """
        plink --bfile $OUTPUT.MERGED.FOUNDERS --freq --out $OUTPUT.MERGED.FOUNDERS.FRQ
        awk '{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > all.remove.snps
        awk '{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}' $OUTPUT.MERGED.FOUNDERS.FRQ.frq > allele.order

        # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
        plink --bfile $OUTPUT.MERGED.FOUNDERS --reference-allele allele.order --exclude all.remove.snps --geno 0.5 --make-bed --out $OUTPUT

        # Calculate allele frequencies
        plink --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ
        rm $SNP_DATA.FOUNDERS.*
        rm $OUTPUT.MERGED.FOUNDERS.*
        rm $OUTPUT.*.FOUNDERS.???
        rm allele.order
        rm all.remove.snps

        """

        TMP_allele_order = os.path.join(INTERMEDIATE_PATH, "allele.order")
        TMP_all_remove_snps = os.path.join(INTERMEDIATE_PATH, "all.remove.snps")

        command = ' '.join([plink, "--bfile", OUTPUT+'.MERGED.FOUNDERS', "--freq", "--out", OUTPUT+'.MERGED.FOUNDERS.FRQ'])
        print(command)
        os.system(command)
        command = ' '.join(["awk", "'{if (NR > 1 && ($5 < 0.0001 || $5 > 0.9999)){print $2}}'", OUTPUT+'.MERGED.FOUNDERS.FRQ.frq', ">", TMP_all_remove_snps])
        print(command)
        os.system(command)
        command = ' '.join(["awk", '\'{if (NR > 1){if (($3 == "A" && $4 == "P") || ($4 == "A" && $3 == "P")){print $2 "\tP"}}}\'', OUTPUT+'.MERGED.FOUNDERS.FRQ.frq', ">", TMP_allele_order])
        print(command)
        os.system(command)

        # QC: Maximum per-SNP missing > 0.5, MAF > 0.1%
        command = ' '.join([plink, "--bfile", OUTPUT+'.MERGED.FOUNDERS', "--reference-allele", TMP_allele_order, "--exclude", TMP_all_remove_snps, "--geno 0.5", "--make-bed", "--out", OUTPUT])
        print(command)
        os.system(command)

        # Calculate allele frequencies
        command = ' '.join([plink, "--bfile", OUTPUT, "--keep-allele-order", "--freq", "--out", OUTPUT+'.FRQ'])
        print(command)
        os.system(command)

        index += 1


        rm_tlist = (SNP_DATA2+'.FOUNDERS.*', OUTPUT+'.MERGED.FOUNDERS.*', OUTPUT+'.*.FOUNDERS.???', TMP_allele_order, TMP_all_remove_snps)

        for i in rm_tlist:
            print(i)
            os.system("rm "+i)



    if PREPARE:

        """
        [Source from Buhm Han.]

        awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
        plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
        awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
        cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped

        echo "[$i] Converting to beagle format.";  @ i++
        linkage2beagle pedigree=$OUTPUT.nopheno.ped data=$OUTPUT.dat beagle=$OUTPUT.bgl standard=true > $OUTPUT.bgl.log


        [Source from Yang.]

        awk '{print $2 " " $4 " " $5 " " $6}' $OUTPUT.bim > $OUTPUT.markers
        plink --bfile $OUTPUT --keep-allele-order --recode --alleleACGT --out $OUTPUT
        plink --bfile $OUTPUT --recode --transpose --out $OUTPUT
        # awk '{print "M " $2}' $OUTPUT.map > $OUTPUT.dat
        # cut -d ' ' -f1-5,7- $OUTPUT.ped > $OUTPUT.nopheno.ped

        echo "[$i] Converting to beagle format.";  @ i++
        beagle2vcf -fnmarker $OUTPUT.markers -fnfam $OUTPUT.fam -fngtype $OUTPUT.tped -fnout $OUTPUT.vcf

        I will make this code block based on source given by Yang. for now.

        """

        """
        [2019. 01. 30.]
        (1) Encode genomic position
        (2) make tped
        (3) Encode alleles
        (4) makebed
        (5) *.markers
        (5.5) recode to ped, map file
        (6) *.dat and *.nopheno.ped => .bgl (by linkage2beagle.jar)
        (7) *.vcf (by beagle2vcf.jar)
        """

        print("[{}] Preparing files for Beagle.".format(index))

        command = ' '.join(["awk", '\'{print $2 " " $4 " " $5 " " $6}\'', OUTPUT+'.bim', ">", OUTPUT+'.markers'])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", OUTPUT, "--keep-allele-order", "--recode", "--alleleACGT", "--out", OUTPUT])
        print(command)
        os.system(command)
        command = ' '.join([plink, "--bfile", OUTPUT, "--recode --transpose --out", OUTPUT])
        #command = ' '.join(["awk", '\'{print "M " $2}\'', OUTPUT+'.map', ">", OUTPUT+'.dat'])
        #print(command)
        os.system(command)
        #command = ' '.join(["cut -d ' ' -f1-5,7-", OUTPUT+'.ped', ">", OUTPUT+'.nopheno.ped'])
        print(command)
        os.system(command)


        index += 1



        print("[{}] Converting to beagle format.".format(index))

        #command = ' '.join([linkage2beagle, "pedigree="+OUTPUT+'.nopheno.ped', "data="+OUTPUT+'.dat', "beagle="+OUTPUT+'.bgl', "standard=true", ">", OUTPUT+'.bgl.log'])
        command = ''.join([beagle2vcf, " -fnmarker ", OUTPUT,".markers -fnfam ", OUTPUT,".fam -fngtype ", OUTPUT,".tped -fnout ", OUTPUT,".vcf"])

        print(command)
        os.system(command)


        index += 1



    if PHASE:

        '''
        beagle unphased=$OUTPUT.bgl nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000 log=$OUTPUT.phasing >> $OUTPUT.bgl.log

        '''
        print("[{}] Phasing reference using Beagle (see progress in $OUTPUT.bgl.log).".format(index))

        #command= ' '.join([beagle, "unphased="+OUTPUT+'.bgl', "nsamples=4 niterations=10 missing=0 verbose=true maxwindow=1000", "log="+OUTPUT+'.phasing', ">>", OUTPUT+'.bgl.log'])

        #new beagle (>v4), assuming 8 threads and 10 interations
        command= ' '.join([beagle, "gt="+OUTPUT+'.vcf', "chrom=6 nthreads=8 niterations=10 lowmem=true", "out="+OUTPUT+'.bgl', ">>", OUTPUT+'.bgl.log'])
        print(command)
        os.system(command)


        index += 1


        #converting back to Beagle v3
        """
        echo "[$i] Converting vcf output to Beagle v3 phased haplotype files."; @ i++

        vcf2phased -gprobs_opt 2 -fnvcf $OUTPUT.bgl.vcf -fnfam $OUTPUT.fam -fnmarker $OUTPUT.markers -fnphased $OUTPUT.bgl.phased

        echo ""
        echo "[$i] Reverting hacked allele identities in *.bgl.vcf file back to those used in reference panel."; @ i++

        revert_alleles -script_opt 2 -fnvcf $OUTPUT.bgl.vcf -fnref $OUTPUT.bim
        """


    if CLEANUP:

        print("[{}] Removing unnecessary files.".format(index))
        '''
        rm $OUTPUT.nopheno.ped
        rm $OUTPUT.bgl.gprobs
        rm $OUTPUT.bgl.r2
        rm $OUTPUT.bgl
        rm $OUTPUT.ped
        rm $OUTPUT.map
        rm $OUTPUT.dat
        rm $OUTPUT.phasing.log
        '''

        rm_tlist = ('.nopheno.ped', '.bgl.gprobs', '.bgl.r2', '.bgl', '.ped', '.map', '.dat')

        for i in rm_tlist:
            print("rm " + OUTPUT + i)
            os.system("rm " + OUTPUT + i)


        index += 1


    print("[{}] Done!".format(index))




    return 0


if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        MakeReference.py

        This script generates a reference panel for HLA imputation

        Usage(1)
        : python3 MakeReference.py 
            --previous-version
            -i ./data/MakeReference_old/HAPMAP_CEU
            -ped ./data/MakeReference_old/HAPMAP_CEU_HLA.ped
            -o ./Trial_HAPMAP_CEU

        Usage(2)
        : python3 MakeReference.py 
            -i ./data/MakeReference/HAPMAP_CEU
            -ped ./data/MakeReference/HAPMAP_CEU_HLA.4field.ped 
            -hg 19 
            -o ./Trial_HAPMAP_CEU
            -dict-AA ./data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt370
            -dict-SNPS ./data/MakeReference/HLA_DICTIONARY_SNPS.hg19.imgt370

        HLA PED file should contain HLA alleles in the following (alphabetical) order:
        HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1

    #################################################################################################
                                     '''),
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-i", help="\nInput Data file prefix(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("-ped", help="\nHLA Type Data(.ped)\n\n", required=True)
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")
    parser.add_argument("-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")
    parser.add_argument("-o", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as original version.\n\n",
                        action='store_true')

    parser.add_argument("-dict-AA", help="\nThe file prefix of HLA Dictionaries of Amino Acids(AA).\n\n")
    # hla_dict.add_argument("-dict-AA-map", help="\nInput HLA Dictionary .map file for AA Information.\n\n", default="Not_given")
    parser.add_argument("-dict-SNPS", help="\nThe file prefix of HLA Dictionaries of SNPs(SNPS).\n\n")
    # hla_dict.add_argument("-dict-SNPS-map", help="\nInput HLA Dictionary .map file for SNPS Information\n\n", default="Not_given")




    ##### <for Test> #####

    # ==========< New version >==========

    # # v370, hg18 (perfectly what Sherman dealt with.)
    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU",
    #                           "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-hg", "18",
    #                           "-o", "./MAKEREFERENCE/MAKEREFERENCE_PYTHON",
    #                           "-dict-AA", "./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-dict-AA-map", "./data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map",
    #                           "-dict-SNPS", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-dict-SNPS-map", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map",
    #                           ])

    # # v3320
    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU",
    #                           "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.ped",
    #                           "-o", "OLD_VERSION_TEST/PREV_VERSION_TEST",
    #                           "-dict-AA", "./data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.txt",
    #                           "-dict-AA-map", "./data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.map",
    #                           "-dict-SNPS", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg19.imgt3320.txt",
    #                           "-dict-SNPS-map", "./data/MakeReference/HLA_DICTIONARY_SNPS.hg19.imgt3320.map",
    #"-m","2000m",
    #                           ])

    # ==========< Perfectly Old version >==========

    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU", "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.ped", "-o", "PREV_VERSION_TEST", "--previous-version"]) # 완전히 구버젼으로 작동하게 하고 싶을때.

    # # Intermediate Path
    # args = parser.parse_args(["-i", "./data/MakeReference/HAPMAP_CEU", "-ped", "./data/MakeReference/HAPMAP_CEU_HLA.ped", "-o", "OLD_VERSION_TEST/PREV_VERSION_TEST", "--previous-version"]) # 완전히 구버젼으로 작동하게 하고 싶을때.



    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)



    ##### Additional Argument processing

    ## Human Genome version

    if args.previous_version:
        args.hg = "18" # Which version is given, override it to "hg 18".


    ## Dictionary

    if not args.previous_version:

        if args.dict_AA and args.dict_SNPS:

            # When all HLA DICTIONARY information is given properly,
            pass

        else:
            # At most one HLA dictionary file prefix is given.
            print(std_ERROR_MAIN_PROCESS_NAME + "The prefix of HLA DICTIONARY files aren't given properly. Please check them all again.")
            print('{"-dict-AA", "-dict-SNPS"}\n')
            sys.exit()



    ## Memory allocation.

    p = re.compile(r'g$')

    if p.search(args.mem):
        args.mem = p.sub(repl="000m", string=args.mem) # Gigabyte to Megabyte to use it in java.
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Given memory for beagle is unappropriate. Please check it again.\n")
        sys.exit()


    # Implementing Main Function.
    MakeReference(_INPUT_DATA=args.i, _hped=args.ped, _OUTPUT=args.o, _dictionary_AA=args.dict_AA,
                  _dictionary_SNPS=args.dict_SNPS, _previous_version=args.previous_version, _hg=args.hg, _mem=args.mem)
