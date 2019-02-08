# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap



########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



def SNP2HLA(_input, _reference_panel, _out,
            _mem = "2000m", _marker_window_size=1000, _tolerated_diff=.15,
            _dependency="./"):



    ### Optional Arguments check.

    if not os.path.exists(_dependency):
        print(std_ERROR_MAIN_PROCESS_NAME + "The path(folder) of depedency('{}') doesn't exist. Please check it again.".format(_dependency))
        sys.exit()

    p_Mb = re.compile(r'\d+m')
    p_Gb = re.compile(r'\d+[gG]')

    if p_Mb.match(_mem):
        pass # No problem.
    elif p_Gb.match(_mem):
        _mem = re.sub(r'[gG]', '000m', _mem)
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Given Java memory value('{}') has bizzare representation. Please check it again.".format(_mem))
        sys.exit()



    ### Check dependencies.

    _plink = os.path.join(_dependency, "plink")         # Plink(v1.9)
    _beagle = os.path.join(_dependency, "beagle.jar")   # Beagle(v4.1)
    _linkage2beagle = os.path.join(_dependency, "linkage2beagle.jar")
    _beagle2linkage = os.path.join(_dependency, "beagle2linkage.jar")
    _merge_table = os.path.join("src/merge_tables.pl")
    _parse_dosage = os.path.join("src/ParseDosage.csh")


    if not os.path.exists(_plink):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare PLINK(v1.90) in 'dependency/' folder.")
        sys.exit()
    if not os.path.exists(_beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare Beagle(v.4.1) in 'dependency/' folder.")
        sys.exit()
    if not os.path.exists(_linkage2beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'linkage2beagle.jar' in 'dependency/' folder.")
        sys.exit()
    if not os.path.exists(_beagle2linkage):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'beagle2linkage.jar' in 'dependency/' folder.")
        sys.exit()
    if not os.path.exists(_merge_table):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'merge_tables.pl' in 'src/' folder.")
        sys.exit()
    if not os.path.exists(_parse_dosage):
        print(std_ERROR_MAIN_PROCESS_NAME + "Please prepare 'ParseDosage.csh' in 'src/' folder.")
        sys.exit()



    ### Intermediate path.

    OUTPUT = _out if not _out.endswith('/') else _out.rstrip('/')
    if bool(os.path.dirname(OUTPUT)):
        INTERMEDIATE_PATH = os.path.dirname(OUTPUT)
        os.makedirs(INTERMEDIATE_PATH, exist_ok=True)
    else:
        # If `os.path.dirname(OUTPUT)` doesn't exist, then it means the output of MakeReference should be genrated in current directory.
        INTERMEDIATE_PATH = "./"


    JAVATMP = _out+".javatmpdir"
    os.system("mkdir -p " + JAVATMP)



    ### Setting commands

    PLINK = ' '.join([_plink, "--silent", "--allow-no-sex"]) # "--noweb" won't be included because it is Plink1.9
    BEAGLE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _beagle])
    LINKAGE2BEAGLE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _linkage2beagle])
    BEAGLE2LINKAGE = ' '.join(["java", "-Djava.io.tmpdir="+JAVATMP, "-Xmx"+_mem, "-jar", _beagle2linkage])



    ### Control Flags
    EXTRACT_MHC = 1
    FLIP = 0
    CONVERT_IN = 0
    IMPUTE = 0
    CONVERT_OUT = 0
    CLEANUP = 0


    print("SNP2HLA: Performing HLA imputation for dataset {}".format(_input))
    print("- Java memory = {}(Mb)".format(_mem))
    print("- Beagle(v4.1) window size = \"{}\" markers".format(_marker_window_size))


    index= 1
    __MHC__ = _out+".MHC"


    if EXTRACT_MHC:

        print("[{}] Extracting SNPs from the MHC.".format(index)); index += 1



    if FLIP:

        print("[{}] Performing SNP quality control.".format(index)); index += 1


    if CONVERT_IN:

        print("[{}] Converting data to beagle format.".format(index)); index += 1


    if IMPUTE:

        print("[{}] Performing HLA imputation.".format(index)); index += 1


    if CONVERT_OUT:

        print("[{}] Converting posterior probabilities to Plink dosage format.".format(index)); index += 1



        print("[{}] Converting imputation genotypes to Plink *.ped format.".format(index)); index += 1



    print("Done\n")

    return 0



if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        < SNP2HLA.py >

        Author: Sherman Jia (xiaomingjia@gmail.com)
                + Small modifications by Buhm Han (buhmhan@broadinstitute.org): 8/7/12
                + Extensive modifications by Phil Stuart (pstuart@umich.edu) on 4/19/16 to allow use of Beagle 4.1:
                  verified to work with 22Feb16, 15Apr16, and 03May16 versions of Beagle.
                + Small modifications by Yang Luo (yangluo@broadinstitute.org): 09/30/16: verfiied working with Bealge 4.1 27Jun16.
                + Recoded to python by Wanson Choi(wschoi.bhlab@gmail.com) : 2019/02/06 
                
                
        DESCRIPTION: This script runs imputation of HLA amino acids and classical alleles using SNP data.        
        
        INPUTS:
        1. Plink dataset (*.bed/bim/fam)
        2. Reference dataset (*.bgl.phased, *.markers in beagle 3.0.4 format; *.fam/.bim/.FRQ.frq in PLINK format)
        
        DEPENDENCIES: (download and place in the same folder as this script)
        1. PLINK (1.9)  (Will not work with older Plink 1.07)
        2. Beagle (4.1) (Need to rename java executable as beagle.jar)
        3. merge_tables.pl (Perl script to merge files indexed by a specific column)
        4. vcf2gprobs.jar (Beagle utility for generating a Beagle v3 genotypes probability file from a Beagle 4.1 vcf file with GT field data)
        5. beagle2vcf and vcf2phased (Utilities by Phil Stuart for vcf <-> Beagle v3 format)
        6. ParseDosage.csh (Converts Beagle posterior probabilities [.gprobs] to dosages in PLINK format [.dosage])
        7. revert_alleles (Utility by Phil Stuart to revert hacked non-ACTG alleles in *.vcf, *.gprobs and *.dosage output files (necessitated by Beagle 4.1) back to allele identities in reference panel)
        8. If genetic_map_file argument is specified, PLINK format genetic map on cM scale (plink.chr6.GRCh36.map, downloaded from http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)
        
        # USAGE: ./SNP2HLA_new.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers/.fam/.bim/.bed/.FRQ.frq) OUTPUT plink max_memory[gb] nthreads niterations genetic_map_file


    #################################################################################################
                                     '''),
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-i", help="\nInput Plink data file prefix(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix\n\n", required=True)


    ##### <for Test> #####


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)



    ##### Additional Argument processing

