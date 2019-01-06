# -*- coding: utf-8 -*-

import os, sys, re
import pandas as pd
import argparse, textwrap


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names2 = [0, 2, 1, 7, 5, 6, 3, 4]  # Read in the order of `HLA_names`, Write in the order of `HLA_names2` (to be compatible with old version of Dictionary).

isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}


def HLAtoSequences(_chped, _dictionary, _type, _out, __use_Pandas=False, __previous_version=False):



    ########## < Argument checking > ##########

    # (1) ped file existence
    if not os.path.isfile(_chped):
        print(std_MAIN_PROCESS_NAME + "Given ped file doen't exist. Please check it again.\n")
        sys.exit()

    # (2) HLA DICTIONARY file
    if not os.path.isfile(_dictionary):
        print(std_MAIN_PROCESS_NAME + "Given dictionary file doen't exist. Please check it againg.\n")
        sys.exit()

    # (3) Chekcing `_type`
    if not (_type == "AA" or _type == "SNPS"):
        print(std_MAIN_PROCESS_NAME + "Given value for argument `_type` has wrong value. Please check it again.\n")
        sys.exit()


    if __use_Pandas:


        ########## < Core Variables > ##########

        HLA_DICTIONARY = pd.DataFrame()
        HLA_DICTIONARY_dict = {}
        ALLELES_SEQ_LENGTH = {}


        ########## < Control Flags > ##########

        LOADING_DICTIONARY = 1
        LOADING_COATED_PED = 1
        BRINGING_SEQUENCES = 1
        EXPORTING_OUTPUT_PED = 1



        if LOADING_DICTIONARY:

            ########## <1. Dictionary Preparation> ##########

            print("\n[1] Loading Dictionary Data.\n")

            if os.path.isfile(_dictionary):
                HLA_DICTIONARY = pd.read_table(_dictionary, sep='\t', header=None, names=["Alleles", "Seqs"], index_col=0)
            else:
                print(std_MAIN_PROCESS_NAME + "Given Dictionary file doesn't exit!\n")
                sys.exit()

            # (2018. 7. 13.) deprecated - going back to use HLA gene captioned way.
            # # Processing HLA gene caption parts(splited by '*')
            # df_temp = pd.DataFrame(HLA_DICTIONARY.loc[:, "Alleles"].apply(lambda x : x.split('*')).tolist(), columns=["HLA", "Alleles"])
            #
            # HLA_DICTIONARY = pd.concat([df_temp, HLA_DICTIONARY.loc[:, "Seqs"]], axis=1).set_index('HLA')
            print(HLA_DICTIONARY.head())


            ### Dividing `HLA_DICTIONARY` in advance.

            # For performance efficiency, `HLA_DICTIONARY` will be divded by HLA gene type in advance.
            HLA_DICTIONARY_dict = {HLA_names[i]: HLA_DICTIONARY.filter(regex= ''.join(["^", HLA_names[i], "\*"]), axis=0) for i in range(0, len(HLA_names))}
            # HLA_DICTIONARY_dict = {HLA_names[i]: HLA_DICTIONARY.loc[HLA_names[i], :].reset_index(drop=True).set_index('Alleles') for i in range(0, len(HLA_names))}

            for i in range(0, len(HLA_names)):
                print("\nSequence Information of %s" % HLA_names[i])
                print(HLA_DICTIONARY_dict[HLA_names[i]].head())

            ALLELES_SEQ_LENGTH = {HLA_names[i] : len(HLA_DICTIONARY_dict[HLA_names[i]].iat[0, 0]) for i in range(0, len(HLA_names))}

            for i in range(0, len(HLA_names)):
                print("\nSequence Length of %s" % HLA_names[i])
                print(ALLELES_SEQ_LENGTH[HLA_names[i]])



        if LOADING_COATED_PED:

            ########## <2. Loading Coated PED(Input PED) file> ##########

            print("\n[2] Loading Input PED file.")

            INPUT_PED = pd.read_table(_chped, sep='\t', header=None, index_col=[0, 1, 2, 3, 4, 5], dtype=str)
            INPUT_PED.columns = pd.Index([name + '_' + str(j + 1) for name in HLA_names for j in range(0, 2)])

            print(INPUT_PED.head())



        if BRINGING_SEQUENCES:

            ########## <3. Bringing Sequences> ##########

            print("\n[3]: Transforming Allele names to Sequences.")

            l_FOR_OUTPUT = []
            l_FOR_OUTPUT_test = []

            for i in range(0, len(HLA_names)):

                curr_dict = HLA_DICTIONARY_dict[HLA_names[i]]

                df_temp = INPUT_PED.filter(regex='_'.join([HLA_names[i], '\d{1}']), axis=1)
                # print(df_temp.head())

                df_temp = df_temp.applymap(lambda x : BringSequence(x, curr_dict) if x != "0" else x)
                # print(df_temp.head())

                l_FOR_OUTPUT_test.append(df_temp)

                # print("\n===============\n")

                # Now, we need to iterate on the number of rows of this DataFrame

                COLUMNS_BY_HLA = []
                # COLUMNS_BY_HLA_test = []

                for j in range(0, len(df_temp)):

                    if df_temp.iat[j, 0] != '0' and df_temp.iat[j, 1] != '0':

                        # (Case1) Most normal case - wehn allele_name is found as pair.
                        # ex) allele1 : A*25:01:01  /  allele2 : A*03:01:01:01

                        # seq1 = df_temp.iat[j, 0] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 0][::-1]
                        # seq2 = df_temp.iat[j, 1] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 1][::-1]

                        # (2018. 3. 9) 다시 여기서 reverse안시키는 걸로 바꿈
                        seq1 = df_temp.iat[j, 0]
                        seq2 = df_temp.iat[j, 1]

                        PAIRED = [value for item in zip(seq1, seq2) for value in item]

                    else:

                        # (Case2) when not found as a pair of alleles, but as a only single allele, => 0 is given
                        # (0, 0) will be compensated as length of `HLA_seq`, Due to here, I need to prepared `len(HLA_seq)` beforehand.

                        seq1 = ['0' for z in range(0, ALLELES_SEQ_LENGTH[HLA_names[i]])]

                        PAIRED = [value for item in zip(seq1, seq1) for value in item]

                    COLUMNS_BY_HLA.append(PAIRED)
                    # COLUMNS_BY_HLA_test.append(''.join(PAIRED))


                l_FOR_OUTPUT.append(pd.DataFrame(COLUMNS_BY_HLA))
                # l_FOR_OUTPUT_test.append(pd.Series(COLUMNS_BY_HLA_test))

                # End of interation. Don't get confused.



        if EXPORTING_OUTPUT_PED:

            ########## <4. Exporting OUTPUT PED file> ##########

            print("\n[4]: Exporting OUTPUT PED file.")

            df_OUTPUT = pd.concat(l_FOR_OUTPUT, axis=1)
            df_OUTPUT.index = INPUT_PED.index
            df_OUTPUT.columns = pd.Index(range(0, df_OUTPUT.shape[1]))

            print(df_OUTPUT.head())



            ### Final Output ped file.
            if _type == 'AA':
                df_OUTPUT.to_csv(_out + '.AA.ped', sep='\t', header=False, index=True)
            elif _type == 'SNPS':
                df_OUTPUT.to_csv(_out + '.SNPS.ped', sep='\t', header=False, index=True)

    else:

        # Processing line by line with python Generators

        ###### < Core Variables > #####

        HLA_DICTIONARY_byHLA = {HLA_names[i]: {} for i in range(0, len(HLA_names))}  # Initialization.
        HLA_SEQ_LENGTH = {HLA_names[i]: -1 for i in range(0, len(HLA_names))}

        HLA_INS_byHLA = {HLA_names[i]: {} for i in range(0, len(HLA_names))}  # Initialization.
        haveInsertion= {HLA_names[i]: -1 for i in range(0, len(HLA_names))}


        if __previous_version:

            ##### < [1] Loading HLA Dictionary > #####

            with open(_dictionary, "r") as f_dictionary:

                count = 0

                for line in f_dictionary:

                    t_line = re.split(r'\s+', line.rstrip('\n'))
                    # ex1) (AA) ['B*58:01:01', 'MLVMAPRTVLLLLSAALALTETWAG...', 'x'] (len == 3)
                    # ex2) (SNPS) ['B*58:01:01', 'CTAGTCCTGCTTCAGGGTCCGGGGCCCG...'] (len == 2)
                    # print(t_line)

                    for i in range(0, len(HLA_names)):

                        if re.match(r'{}:'.format(HLA_names[i]), t_line[0]):

                            HLA_DICTIONARY_byHLA[HLA_names[i]][t_line[0]] = t_line[1] # Sequence information.

                            if HLA_SEQ_LENGTH[HLA_names[i]] == -1:
                                HLA_SEQ_LENGTH[HLA_names[i]] = len(t_line[1])


                            if _type == "AA":

                                HLA_INS_byHLA[HLA_names[i]][t_line[0]] = t_line[2] # Insertion information.

                                if haveInsertion[HLA_names[i]] == -1:
                                    haveInsertion[HLA_names[i]] = bool(t_line[2])

                            break  # If a line of given dictionary belongs to either HLA, then checking whether it belongs to other HLAs is useless.

                    count += 1
                    # if count > 5 : break

            # # Result check
            # for i in range(0, len(HLA_names)):
            #     print("\n{} :\n".format(HLA_names[i]))
            #     for k, v in HLA_DICTIONARY_byHLA[HLA_names[i]].items():
            #         print("{} : {}".format(k, v))
            #
            # for k, v in HLA_SEQ_LENGTH.items():
            #     print("The length of HLA-{} : {}".format(k, v))
            #
            # print("Insertion check : {}".format(haveInsertion))


            ##### < [2] Transforming each HLA alleles to corresponding sequences > #####

            with open(_out + ".{}.ped".format(_type), 'w') as f_output:
                f_output.writelines(GenerateLines(_chped, _type, HLA_DICTIONARY_byHLA, HLA_SEQ_LENGTH, HLA_INS_byHLA, haveInsertion))


        else:
            print(std_MAIN_PROCESS_NAME + "This function will be completed a little bit later.")


def BringSequence(_single_allele, _dict):

    try:
        Seq = _dict.loc[_single_allele, "Seqs"]
    except KeyError:
        Seq = "0"

    return Seq


def GenerateLines(_chped, _type, _dict_seq, _seq_length, _dict_ins, _haveIns):

    with open(_chped, "r") as f:

        for line in f:
            t_line = re.split(r'\s+', line.rstrip('\n'))
            # print(t_line)

            """
            [0,1,2,3,4,5] := ped file information
            [6,7] := HLA-A,
            [8,9] := HLA-B,
            ...,
            [20, 21] := HLA-DRB1
            """

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_part__ = '\t'.join([
                BringSequence2(t_line[(2 * i + 6)], t_line[(2 * i + 7)], _type, HLA_names[i],
                               _dict_seq[HLA_names[i]], _seq_length[HLA_names[i]],
                               _dict_ins[HLA_names[i]], _haveIns) for i in HLA_names2
            ])

            # mem_p2 = process.memory_info().rss / 1024 ** 2
            # print("{}(Mb)".format(mem_p2 - mem_p1))

            yield '\t'.join([__ped_info__, __genomic_part__]) + "\n"



def BringSequence2(_HLA_allele1, _HLA_allele2, _type, _hla,
                   _dict_seq, _seq_length,
                   _dict_ins, _haveIns):

    # # checking `_dict`
    # count = 0
    #
    # for k,v in _dict.items():
    #     print("{} : {}".format(k, v))
    #
    #     count += 1
    #     if count > 10:
    #         break

    # print("al1 : {}\nal2 : {}".format(_HLA_allele1, _HLA_allele2))

    # Finding the corresponding sequence and insertion marker of `_HLA_allele1`.
    try:
        Seq1 = _dict_seq[_HLA_allele1]

        if _type == "AA" and isREVERSE[_hla]:
            Seq1 = Seq1[::-1]

    except KeyError:
        Seq1 = "-1" # fail


    # Same job for `_HLA_allele2`.
    try:
        Seq2 = _dict_seq[_HLA_allele2]

        if _type == "AA" and isREVERSE[_hla]:
            Seq2 = Seq2[::-1]

    except KeyError:
        Seq2 = "-1" # fail


    # print("Corresponding Seqs : \n{}\n{}".format(Seq1, Seq2))

    ### Main sequence information processing
    if Seq1 != "-1" and Seq2 != "-1":

        # Only when both HLA alleles can get the corresponding HLA sequence information.

        l_temp = []

        for i in range(0, len(Seq1)):
            l_temp.append(Seq1[i])
            l_temp.append(Seq2[i])

        # Reversing
        __return__ = '\t'.join(l_temp)

    else:
        __return__ = '\t'.join(["0" for z in range(0, 2 * _seq_length)])



    ### Insertion part processing.
    if _type == "AA" and _haveIns[_hla]:

        # One important assumption : "Every insertion part is just one, single character.(ex 'x', 'P', 'A')"
        try:
            ins1 = _dict_ins[_HLA_allele1]
        except KeyError:
            ins1 = ""

        try:
            ins2 = _dict_ins[_HLA_allele2]
        except KeyError:
            ins2 = ""

        if bool(ins1) and bool(ins2):
            __insertions__ = '\t'.join([ins1, ins2])
        else:
            __insertions__ = '\t'.join(["0", "0"])


        __return__ = '\t'.join([__return__, __insertions__])




    return __return__




if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

        Original Author: Sherman Jia, 2012

        HLAtoSequences.py
        - This script Converts HLA alleles (in .ped file format) to amino acid or DNA sequences

        Input file should contain: FID, IID, pID, mID, sex, pheno, HLA-A (2), B (2), C (2),
        DPA1 (2), DPB1 (2), DQA1 (2), DQB1 (2), DRB1 (2) ... Broad Order

    #########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"
    # parser._optionals.description = "- Necessary main options.\n"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-dict", help="\nHLA dictonary file name(ex. 'HLA_DICTIONARY_AA.txt')\n\n", required=True)
    parser.add_argument("-type", help="\nAA(for Amino Acid) or SNP(for SNPs)\n\n", choices=["AA", "SNPS"], required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as original version.\n\n",
                        action='store_true')




    ##### <for Test> #####

    # (2018.2.9)
    # args = parser.parse_args(["./HAPMAP_CEU_HLA_switched.ped", "./HLA_DICTIONARY_AA_hg19.txt", "AA", "--OUTPUT", "BROAD_ORDER"])


    # (2018.2.26)
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./HLA_DICTIONARY_AA.hg19.imgt370.txt", "AA", "--OUTPUT", "TEST_0329_HLAtoSeq"])

    # (2018.3.8)
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./HLA_DICTIONARY_SNPS.hg19.imgt370.txt", "SNPS", "--OUTPUT", "BROAD_ORDER"])

    # (2018. 7. 12.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370"])



    ## (2018. 01. 06.) Introducinig compatibility to work with old version of dictionary

    # # AA
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference_old/HLA_DICTIONARY_AA.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/HAPMAP_CEU.old.enCODED",
    #                           "--previous-version"])

    # # SNPS
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
    #                           "-dict", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference_old/HLA_DICTIONARY_SNPS.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/HAPMAP_CEU.old.enCODED",
    #                           "--previous-version"])


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)


    HLAtoSequences(args.ped, args.dict, args.type, args.o, __previous_version=args.previous_version)
