# -*- coding: utf-8 -*-

import os, re
import pandas as pd
import argparse, textwrap


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


def encodeVariants(_ped, _map, _out, __use_pandas=False, __asCapital=True, __addDummyMarker=False):


    if __use_pandas:

        ########## < 1. Loading MAP file > ##########

        print(std_MAIN_PROCESS_NAME + "[1] Loading MAP file.\n")

        # Processing map file
        df_mapfile = pd.read_table(_map, sep='\t', header=None, dtype=str)
        print(df_mapfile.head())


        ########## < 2. Allele overlapping > ##########

        print(std_MAIN_PROCESS_NAME + "[2] Allele Overlapping.\n")

        df_pedfile = pd.read_table(_ped, sep='\t', header=None, dtype=str)
        print(df_pedfile.head())

        """
        # Before using "apply()" function
        alleles = [[] for i in range(0, int((len(df_pedfile.columns)-6)/2))] # 2396
    
        for i in range(0, len(df_pedfile.index)):
    
            for j in range(6, len(df_pedfile.columns), 2):
    
                SNP1 = df_pedfile.iat[i, j]
                SNP2 = df_pedfile.iat[i, j+1]
    
                idx_alleles = int((j - 6)/2)
    
                if SNP1 != "0":
                    if SNP1 not in alleles[idx_alleles]:
                        # alleles[idx_alleles] += [SNP1]
                        alleles[idx_alleles].extend((SNP1,))
    
                if SNP2 != "0":
                    if SNP2 not in alleles[idx_alleles]:
                        # alleles[idx_alleles] += [SNP2]
                        alleles[idx_alleles].extend((SNP2,))
    
                # print("%s %s" % (SNP1, SNP2))
        """

        flattened = df_pedfile.iloc[:, 6:].apply(set, axis=0)
        alleles = [tuple(flattened.iat[i].union(flattened.iat[i+1]).difference({0, "0"})) for i in range(0, len(flattened), 2)]

        print(alleles)
        # Now, `alleles`(list of list) will be used to generate new .ped file.
        print("Making Alleles list is done!")


        ########## < 3. Making new .ped file > ##########

        print(std_MAIN_PROCESS_NAME + "[3] Making new .ped file.\n")

        # each lines will be converted to lists,
        # and these lists will be gathered to make list of lists. Finally this list of lists will be DataFrame.

        to_DataFrame = []

        # Iterating over each lines of .ped files
        for N_pedline in range(0, len(df_pedfile.index)):

            curr_line = list(df_pedfile.iloc[N_pedline, :])

            ID = curr_line[0:6]
            Seq = []

            for i in range(6, len(curr_line), 2):

                # (SNP1, SNP2) : (6,7) -> (8,9) -> (10, 11) -> ...

                SNP1 = curr_line[i]
                SNP2 = curr_line[i+1]

                idx_alleles = int((i - 6)/2)
                # idx_Seq = i-6

                # alleles[idx_alleles] := the number of sort of alleles on each positions.
                if len(alleles[idx_alleles]) > 2:

                    for j in range(0, len(alleles[idx_alleles])):

                        if SNP1 == "0" or SNP2 == "0":
                            Seq.append("0"); Seq.append("0")

                        else:

                            if alleles[idx_alleles][j] == SNP1:
                                Seq.append("P")
                            else:
                                Seq.append("A")

                            if alleles[idx_alleles][j] == SNP2:
                                Seq.append("P")
                            else:
                                Seq.append("A")


                    if len(alleles[idx_alleles]) > 3:

                        j_end = 1 if len(alleles[idx_alleles]) == 4 else len(alleles[idx_alleles])

                        for j in range(0, j_end):

                            for k in range(j+1, len(alleles[idx_alleles])):

                                if SNP1 == "0" or SNP2 == "0":
                                    Seq.append("0"); Seq.append("0")

                                else:
                                    if alleles[idx_alleles][j] == SNP1 or alleles[idx_alleles][k] == SNP1:
                                        Seq.append("P")
                                    else:
                                        Seq.append("A")

                                    if alleles[idx_alleles][j] == SNP2 or alleles[idx_alleles][k] == SNP2:
                                        Seq.append("P")
                                    else:
                                        Seq.append("A")


                        if len(alleles[idx_alleles]) > 5:

                            j_end = 1 if len(alleles[idx_alleles]) == 6 else len(alleles[idx_alleles])

                            for j in range(0, j_end):
                                for k in range(j+1, len(alleles[idx_alleles])):
                                    for l in range(k+1, len(alleles[idx_alleles])):

                                        if SNP1 == "0" or SNP2 == "0":
                                            Seq.append("0"); Seq.append("0")

                                        else:
                                            if alleles[idx_alleles][j] == SNP1 or alleles[idx_alleles][k] == SNP1 or alleles[idx_alleles][l] == SNP1:
                                                Seq.append("P")
                                            else:
                                                Seq.append("A")

                                            if alleles[idx_alleles][j] == SNP2 or alleles[idx_alleles][k] == SNP2 or alleles[idx_alleles][l] == SNP2:
                                                Seq.append("P")
                                            else:
                                                Seq.append("A")



                else:
                    # Most of Cases have length less than equal 2, they will fall into this if-else block.
                    if SNP1 == "0" or SNP2 == "0":
                        Seq.append("0"); Seq.append("0")
                    else:
                        Seq.append(curr_line[i]); Seq.append(curr_line[i+1])


            to_DataFrame.append(ID+Seq)

            # End-point of one line of .ped file.

        df_output_pedfile = pd.DataFrame(to_DataFrame)
        df_output_pedfile.to_csv(_out + '.ped', sep='\t', header=False, index=False)

        print("Making .ped file is done")


        ########## < 4. Making new .map file > ##########

        print(std_MAIN_PROCESS_NAME + "[4] Making new .map file.\n")

        to_DataFrame = []

        # print(len(alleles))
        # print(len(df_mapfile.index))
        # len(alleles) == len(df_mapfile.index)

        for i in range(0, len(alleles)):

            # print(alleles[i])
            # print(list(df_mapfile.iloc[i, :]))

            curr_line = tuple(df_mapfile.iloc[i, :])
            # [0]: 6, [1]: 'AA_A_9_30126516', [2]: 0, [3]: 30126516

            if len(alleles[i]) > 2:

                # multi_alleles_1,2,3 => These are all list of lists

                multi_alleles_1 = [(curr_line[0], curr_line[1]+'_'+alleles[i][j], curr_line[2], curr_line[3]) for j in range(0, len(alleles[i]))]
                # to_DataFrame += multi_alleles_1
                to_DataFrame.extend(multi_alleles_1)

                if len(alleles[i]) > 3:

                    j_end = 1 if len(alleles[i]) == 4 else len(alleles[i])

                    # for j in range(0, j_end):
                    #     for k in range(j+1, len(alleles[i])):
                    #         print([curr_line[0], curr_line[1]+ '_' + alleles[i][j]+alleles[i][k], curr_line[2], curr_line[3]])

                    multi_alleles_2 = [(curr_line[0], curr_line[1]+ '_' + alleles[i][j]+alleles[i][k], curr_line[2], curr_line[3]) for j in range(0, j_end) for k in range(j+1, len(alleles[i]))]
                    # to_DataFrame += multi_alleles_2
                    to_DataFrame.extend(multi_alleles_2)

                    if len(alleles[i]) > 5:

                        j_end = 1 if len(alleles[i]) == 6 else len(alleles[i])

                        # for j in range(0, j_end):
                        #     for k in range(j+1, len(alleles[i])):
                        #         for l in range(k+1, len(alleles[i])):

                        multi_alleles_3 = [(curr_line[0], curr_line[1]+ '_' + alleles[i][j]+alleles[i][k]+alleles[i][l], curr_line[2], curr_line[3]) for j in range(0, j_end) for k in range(j+1, len(alleles[i])) for l in range(k+1, len(alleles[i]))]
                        # to_DataFrame += multi_alleles_3
                        to_DataFrame.extend(multi_alleles_3)

            # Job to divide multi-allele is done.


            else:
                to_DataFrame.extend([curr_line])

        df_output_mapfile = pd.DataFrame(to_DataFrame)
        df_output_mapfile.to_csv(_out + '.map', sep='\t', header=False, index=False)

        # pd.Series(alleles).to_csv('alleles.txt', sep='\t', header=False)

        # print(df_output_mapfile.head())

        print("Making new .map file is done!")


    else:

        # Processing line by line with python Generators to save memory

        ########## < Control Flags > ##########

        _1_ALLELE_OVERLAPPING = 1
        _2_MAKING_NEW_PEDFILE = 1
        _3_MAKING_NEW_MAPFILE = 1
        _4_MAKING_ALLELELIST = 1


        if _1_ALLELE_OVERLAPPING:

            ########## < [1] Allele overlapping > ##########

            # Acquiring column number
            n_loci = 0
            n_row = 0

            with open(_ped, 'r') as f:

                for line in f:

                    t_line = re.split(r'\s+', line.rstrip('\n'))
                    # print(t_line[6:])

                    if n_row == 0:

                        genomic_info = t_line[6:]
                        n_loci = int(len(genomic_info)/2) # == len(l_factors)

                        # Initializing the list which contains factors which appear in each locus.
                        l_factors = [[] for i in range(0, n_loci)] # Initialization
                        # print(l_factors)


                    for i in range(0, n_loci):

                        idx1 = 2*i + 6 # index for `t_line`
                        idx2 = idx1 + 1

                        ##### Allele overlapping
                        if t_line[idx1] != "0" and (t_line[idx1] not in l_factors[i]):
                            l_factors[i].append(t_line[idx1])
                        if t_line[idx2] != "0" and (t_line[idx2] not in l_factors[i]):
                            l_factors[i].append(t_line[idx2])


                    n_row += 1
                    # if n_row > 5 : break


            # Sorting elements of each lists.
            for i in range(0, len(l_factors)):
                l_factors[i].sort()

            ### --- `l_factors` done.




        if _2_MAKING_NEW_PEDFILE:

            ########## < [2] Making new .ped file > ##########

            with open(_out + ".ped", 'w') as f_NewPed:
                f_NewPed.writelines(MakeNewPed(_ped, l_factors, __asCapital, __addDummyMarker))



        if _3_MAKING_NEW_MAPFILE:

            ########## < [3] Making new .map file > ##########

            with open(_out + ".map", 'w') as f_NewMap:
                f_NewMap.writelines(MakeNewMap(_map, l_factors, __addDummyMarker))



        if _4_MAKING_ALLELELIST:

            ########## < [4] Making *.allelelist file > ##########

            with open(_out + ".allelelist", 'w') as f_allelelist:
                f_allelelist.writelines(MakeAlleleList(_map, l_factors))




def divideToBinaryMarkers(_SNP1, _SNP2, _factors, __asCapital=True):

    _present_ = "P" if __asCapital else "p"
    _absent_ = "A" if __asCapital else "a"

    Seq = []

    if len(_factors) > 2:

        for j in range(0, len(_factors)):

            if _SNP1 == "0" or _SNP2 == "0":
                Seq.append("0"); Seq.append("0")

            else:

                if _factors[j] == _SNP1:
                    Seq.append(_present_)
                else:
                    Seq.append(_absent_)

                if _factors[j] == _SNP2:
                    Seq.append(_present_)
                else:
                    Seq.append(_absent_)

        if len(_factors) > 3:

            j_end = 1 if len(_factors) == 4 else len(_factors)

            for j in range(0, j_end):

                for k in range(j + 1, len(_factors)):

                    if _SNP1 == "0" or _SNP2 == "0":
                        Seq.append("0"); Seq.append("0")

                    else:
                        if _factors[j] == _SNP1 or _factors[k] == _SNP1:
                            Seq.append(_present_)
                        else:
                            Seq.append(_absent_)

                        if _factors[j] == _SNP2 or _factors[k] == _SNP2:
                            Seq.append(_present_)
                        else:
                            Seq.append(_absent_)

            if len(_factors) > 5:

                j_end = 1 if len(_factors) == 6 else len(_factors)

                for j in range(0, j_end):
                    for k in range(j + 1, len(_factors)):
                        for l in range(k + 1, len(_factors)):

                            if _SNP1 == "0" or _SNP2 == "0":
                                Seq.append("0"); Seq.append("0")

                            else:
                                if _factors[j] == _SNP1 or _factors[k] == _SNP1 or _factors[l] == _SNP1:
                                    Seq.append(_present_)
                                else:
                                    Seq.append(_absent_)

                                if _factors[j] == _SNP2 or _factors[k] == _SNP2 or _factors[l] == _SNP2:
                                    Seq.append(_present_)
                                else:
                                    Seq.append(_absent_)



    else:
        # Most of Cases have length less than equal 2, they will fall into this if-else block.
        if _SNP1 == "0" or _SNP2 == "0":
            Seq.append("0"); Seq.append("0")
        else:
            Seq.append(_SNP1); Seq.append(_SNP2)

    return '\t'.join(Seq)



def MakeNewPed(_p_ped, _l_factors, __asCapital=True, __addDummyMarker=False):

    count = 0

    with open(_p_ped, 'r') as f:
        for line in f:
            t_line = re.split(r'\s+', line.rstrip('\n'))

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_info__ = '\t'.join([
                divideToBinaryMarkers(t_line[2 * i + 6], t_line[2 * i + 7], _l_factors[i], __asCapital) for i in range(0, len(_l_factors))
            ])

            __return__ = '\t'.join([__ped_info__, __genomic_info__])


            if __addDummyMarker:
                # add Dummy Markers.
                dummy_markers = '\t'.join(['d', 'D'] if bool(count % 2) else ['D', 'd'])
                __return__ = '\t'.join([__return__, dummy_markers])


            yield __return__ + "\n"




def MakeNewMap(_p_map, _l_factors, __addDummyMarker=False):

    count = 0

    with open(_p_map, 'r') as f_map:

        for l in f_map:

            idx = count # index for `l_factors`

            if len(_l_factors[idx]) > 2:

                t_line = re.split(r'\s+', l.rstrip('\n'))

                for j in range(0, len(_l_factors[idx])):
                    yield '\t'.join([t_line[0], t_line[1] + '_' + _l_factors[idx][j], t_line[2], t_line[3]]) + "\n"


                if len(_l_factors[idx]) > 3:

                    j_end = 1 if len(_l_factors[idx]) == 4 else len(_l_factors[idx])

                    for j in range(0, j_end):
                        for k in range(j+1, len(_l_factors[idx])):
                            yield '\t'.join([t_line[0], t_line[1]+ '_' + _l_factors[idx][j]+_l_factors[idx][k], t_line[2], t_line[3]]) + "\n"


                    if len(_l_factors[idx]) > 5:

                        j_end = 1 if len(_l_factors[idx]) == 6 else len(_l_factors[idx])

                        for j in range(0, j_end):
                            for k in range(j+1, len(_l_factors[idx])):
                                for l in range(k+1, len(_l_factors[idx])):
                                    yield '\t'.join([t_line[0], t_line[1]+ '_' + _l_factors[idx][j]+_l_factors[idx][k]+_l_factors[idx][l], t_line[2], t_line[3]]) + "\n"

            else:
                yield l # "\n" is included.

            count += 1


        if __addDummyMarker:
            # Adding Dummy Marker
            yield '\t'.join(["6", "DummyMarker", "0", "33999999"]) + "\n"



def MakeAlleleList(_p_map, _l_factors):

    count = 0

    with open(_p_map, 'r') as f_map:

        for l in f_map:

            idx = count

            t_line = re.split(r'\s+', l.rstrip('\n'))

            locus_label = t_line[1]
            alleleset = _l_factors[idx]

            yield '\t'.join([locus_label, str(alleleset)]) + "\n"

            count += 1




if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        Original Author : Sherman Jia, 2012
     
        encodeVariants.py
     
        - This script encodes multi-allelic PLINK markers (amino acids and SNPs) into bi-allelic 
            markers 
        - Input files include a normal PLINK .ped and .map file (where the .ped file contains 
            multi-allelic positions).

    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-map", help="\nMap file for given \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)


    ##### <for Test> #####

    # (2018. 2. 26)
    # args = parser.parse_args(["TEST_v2.AA.ped", "TEST_v2.AA.map", "TEST_v2"])

    # (2018. 7. 13.)

    # # AA
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_AA.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.AA.enCODED"])

    # # SNPS
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HLA_DICTIONARY_SNPS.hg18.imgt370.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370.SNPS.enCODED"])

    # (2019. 1. 6.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/_Case_HAPMAP_CEU.old.AA.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/HLA_DICTIONARY_AA.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_2_encodeVariants/_Case_HAPMAP_CEU.AA.enCODED"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/_Case_HAPMAP_CEU.old.SNPS.ped",
    #                           "-map", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_1_HLAtoSequences/HLA_DICTIONARY_SNPS.map",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_2_encodeVariants/_Case_HAPMAP_CEU.SNPS.enCODED"])


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)


    encodeVariants(args.ped, args.map, args.o)