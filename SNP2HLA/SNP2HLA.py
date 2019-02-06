# -*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap



########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



def SNP2HLA():


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

    parser.add_argument("-i", help="\nInput Data file prefix(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix\n\n", required=True)


    ##### <for Test> #####


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)



    ##### Additional Argument processing

