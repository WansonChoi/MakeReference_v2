# -*- coding: utf-8 -*-

# (2017/11/27) recoded by Wanson Choi
import os, re
import pandas as pd
import argparse, textwrap
from collections import OrderedDict


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names2 = [0, 2, 1, 7, 5, 6, 3, 4]  # Read in the order of `HLA_names`, Write in the order of `HLA_names2` (to be compatible with old version of Dictionary).

# (2018. 9. 25.) Replaced by lift-over values.
genepos_hg = {"18": {"A": 30018226, "C": 31344505, "B": 31429628, "DRB1": 32654525, "DQA1": 32713161, "DQB1": 32735219,
                     "DPA1": 33140324, "DPB1": 33151681},
              "19": {"A": 29910247, "C": 31236526, "B": 31321649, "DRB1": 32546547, "DQA1": 32605183, "DQB1": 32627241,
                     "DPA1": 33032346, "DPB1": 33043703},
              "38": {"A": 29942470, "C": 31268749, "B": 31353872, "DRB1": 32578770, "DQA1": 32637406, "DQB1": 32659464,
                     "DPA1": 33064569, "DPB1": 33075926}}




def encodeHLA(_CHPED, _OUTPUT, _hg="18", __use_Pandas=False, __asCapital=True, __addDummyMarker=False,
              __previous_version=False):





    if __use_Pandas:


        ########## <Core Variables> ##########

        # 주어진 ped파일에 HLA column 별로 나타나는 allele name들의 집합
        ALLELE_TABLES = OrderedDict()
        ALLELE_TABLES_1field = OrderedDict()

        # map파일의 Lable로 준비시킬 항목들 집합
        ALL_ALLELES = []  # 결국 ALL_ALLELES == map_LABELS 라고 생각해도 괜춘.
        dict_ALL_ALLELES = {}  # 추후에 search할때 빨리 하려고

        df_OUTPUT_map = pd.DataFrame()
        df_OUTPUT_ped = pd.DataFrame()


        ########## < Control Flags > ##########


        LOADING_PEDFILE = 1
        MAKING_ALLELE_TABLE = 1
        MAKING_OUTPUT_MAP = 1
        MAKING_OUTPUT_PED = 1



        if LOADING_PEDFILE:

            ########## < 1. Loading Input PED file > ##########

            print(std_MAIN_PROCESS_NAME + "[1] Loading Input PED file.\n")

            INPUT_PED = pd.read_table(_CHPED, sep='\t', header=None, dtype=str,
                                      names = ['Fam_ID', 'Sample_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phe'] + [''.join([HLA_names[i], '_', str(j)]) for i in range(0, len(HLA_names)) for j in range(1,3)],
                                      index_col=[0, 1, 2, 3, 4, 5]
                                      )

            print(INPUT_PED.head())



        if MAKING_ALLELE_TABLE:

            ########## < 2. Making Allele Table > ##########

            print(std_MAIN_PROCESS_NAME + "[2] Making Allele Table.\n")

            """
            for i in range(0, len(INPUT_PED.index)):
    
                line = tuple(INPUT_PED.iloc[i, :])
                # 0:6(6컬럼) => .ped sample header
                # 6: => HLA allele pair
    
                alleles["HLA_A_"+line[6]] = genepos["HLA_A"];       alleles["HLA_A_"+line[7]] = genepos["HLA_A"];
                alleles["HLA_A_"+line[6][0:2]] = genepos["HLA_A"];  alleles["HLA_A_"+line[7][0:2]] = genepos["HLA_A"];
    
                alleles["HLA_B_"+line[8]] = genepos["HLA_B"];       alleles["HLA_B_"+line[9]] = genepos["HLA_B"];
                alleles["HLA_B_"+line[8][0:2]] = genepos["HLA_B"];  alleles["HLA_B_"+line[9][0:2]] = genepos["HLA_B"];
    
                alleles["HLA_C_"+line[10]] = genepos["HLA_C"];       alleles["HLA_C_"+line[11]] = genepos["HLA_C"];
                alleles["HLA_C_"+line[10][0:2]] = genepos["HLA_C"];  alleles["HLA_C_"+line[11][0:2]] = genepos["HLA_C"];
    
                alleles["HLA_DPA1_"+line[12]] = genepos["HLA_DPA1"];       alleles["HLA_DPA1_"+line[13]] = genepos["HLA_DPA1"];
                alleles["HLA_DPA1_"+line[12][0:2]] = genepos["HLA_DPA1"];  alleles["HLA_DPA1_"+line[13][0:2]] = genepos["HLA_DPA1"];
    
                alleles["HLA_DPB1_"+line[14]] = genepos["HLA_DPB1"];       alleles["HLA_DPB1_"+line[15]] = genepos["HLA_DPB1"];
                alleles["HLA_DPB1_"+line[14][0:2]] = genepos["HLA_DPB1"];  alleles["HLA_DPB1_"+line[15][0:2]] = genepos["HLA_DPB1"];
    
                alleles["HLA_DQA1_"+line[16]] = genepos["HLA_DQA1"];       alleles["HLA_DQA1_"+line[17]] = genepos["HLA_DQA1"];
                alleles["HLA_DQA1_"+line[16][0:2]] = genepos["HLA_DQA1"];  alleles["HLA_DQA1_"+line[17][0:2]] = genepos["HLA_DQA1"];
    
                alleles["HLA_DQB1_"+line[18]] = genepos["HLA_DQB1"];       alleles["HLA_DQB1_"+line[19]] = genepos["HLA_DQB1"];
                alleles["HLA_DQB1_"+line[18][0:2]] = genepos["HLA_DQB1"];  alleles["HLA_DQB1_"+line[19][0:2]] = genepos["HLA_DQB1"];
    
                alleles["HLA_DRB1_"+line[20]] = genepos["HLA_DRB1"];       alleles["HLA_DRB1_"+line[21]] = genepos["HLA_DRB1"];
                alleles["HLA_DRB1_"+line[20][0:2]] = genepos["HLA_DRB1"];  alleles["HLA_DRB1_"+line[21][0:2]] = genepos["HLA_DRB1"];
    
    
            for k,v in alleles.items():
                print("key : {0} / value : {1}".format(k, v))
                
            """

            # In the past, `ALLELE_TABLES` was created based on dictionary data structure, but now it is created by "apply()" function with "set()" function.

            for i in range(0, len(HLA_names)):
                # for i in range(0, 1):

                temp = INPUT_PED.filter(regex=HLA_names[i]+'_\d', axis=1).apply(set, axis=0).apply(lambda x : x.difference({0, "0"}))

                # print(temp)

                # Column 1 and 2
                set_al1 = temp.iat[0] # ex)     A_1 := {A*32:01:01, A*23:01:01, A*26:01:01, A*02:01:0...
                set_al2 = temp.iat[1] # ex)     A_2 := {A*24:02:01:01, A*01:01:01:01, A*02:05:01, A*3...

                # sr_Unioned_Set = pd.Series(list(set_al1.union(set_al2))).sort_values() # sorting

                l_Unioned_Set = list(set_al1.union(set_al2))
                l_Unioned_Set.sort() # sorting

                # l_Unioned_Set = pd.Series(l_Unioned_Set)
                print(l_Unioned_Set)

                ALLELE_TABLES[HLA_names[i]] = l_Unioned_Set


                ##### Dealing with 1-field #####

                if len(l_Unioned_Set) > 0:

                    ### The case where the union of set of each two column has at least 1 element.

                    sr_temp_1field = pd.Series(l_Unioned_Set).apply(lambda x : re.match(pattern='\*'.join([HLA_names[i], '\d{2,3}']), string=x).group()).unique() # Going through "unique()" function.

                    # print("\nsr_temp_1field\n")
                    # print(sr_temp_1field)

                    ALLELE_TABLES_1field[HLA_names[i]] = sr_temp_1field.tolist()

                else:
                    ### 집합 원소의 개수가 0 일때,(아예 input_ped에서 allele_name이 '0'으로 주어져서 없는 경우)
                    ALLELE_TABLES_1field[HLA_names[i]] = l_Unioned_Set


            # for k, v in ALLELE_TABLES.items():
            #
            #     print("\n===============\n")
            #     print("{0} : \n{1}".format(k, v))



        if MAKING_OUTPUT_MAP:

            ########## < 3. Making OUTPUT .map file > ##########

            print(std_MAIN_PROCESS_NAME + "[3] Making OUTPUT .map file.\n")

            """        
            to_df_OUTPUT_map = []
    
    
            for name in HLA_names:
                for k in sorted_keys:
                    temp = k.split('_')
                    al = temp[2]
    
                    if (temp[1] == name) and (al != "NA") and (al != "") and (al != "0") and (al != "0 0"):
                        # to_df_OUTPUT_map += [ ["6", k, "0", genepos['_'.join([temp[0], temp[1]])]] ]
                        to_df_OUTPUT_map.extend([("6", k, "0", genepos['_'.join([temp[0], temp[1]])])])
    
            df_OUTPUT_map = pd.DataFrame(to_df_OUTPUT_map)
            df_OUTPUT_map.to_csv(_OUTPUT + '.HLA.map', sep='\t', header=False, index=False)
                    
            """


            """ 
            map_Label 만들기 
            
            (1) map_LABELS
            (2) map_CHR
            (3) map_GENETIC_DISTANCE
            (4) map_POS
            
            """

            ##### Making Label for *.map file. #####

            # ALL_ALLELES = pd.concat([ALLELE_TABLES[HLA_names[i]].append(ALLELE_TABLES_1field[HLA_names[i]]) for i in range(0, len(HLA_names))]).sort_values()

            ALL_ALLELES = ALLELE_TABLES[HLA_names[0]]
            ALL_ALLELES.extend(ALLELE_TABLES_1field[HLA_names[0]])
            # ALL_ALLELES = [ALLELE_TABLES[HLA_names[i]].append(ALLELE_TABLES_1field[HLA_names[i]]) for i in range(0, len(HLA_names))]

            for i in range(1, len(HLA_names)):
                ALL_ALLELES.extend(ALLELE_TABLES[HLA_names[i]])
                ALL_ALLELES.extend(ALLELE_TABLES_1field[HLA_names[i]])

            ALL_ALLELES.sort()

            print("\nALL_ALLELES\n")
            print(ALL_ALLELES)


            ### HLA_index(Mining HLA gene name)
            sr_HLA = pd.Series(ALL_ALLELES).apply(lambda x : re.search(pattern='\w+\*', string=x).group().rstrip('*'))
            # print(sr_HLA)

            ### map_LABELS & map_POS ###
            map_LABELS = pd.Series(['HLA_' + ALL_ALLELES[i] for i in range(0, len(ALL_ALLELES))])
            # print("\n`map_LABELS`\n")
            # print(map_LABELS)

            map_POS = [genepos_hg[_hg][sr_HLA.iat[i]] for i in range(0, len(sr_HLA))]
            # print(map_POS)


            # 추후 ped파일에서 search할때 편하게 HLA_name별로 나눠놓기
            dict_ALL_ALLELES = {HLA_names[i] : [ALL_ALLELES[j] for j in range(0, len(ALL_ALLELES)) if (HLA_names[i]+'*' in ALL_ALLELES[j])] for i in range(0, len(HLA_names))}

            # print("\nsegmented `ALL_ALLELES`\n")
            # for k,v in dict_ALL_ALLELES.items():
            #     print("\n============\n")
            #     print("{0} : \n{1}".format(k, v))


            ##### map_Label을 제외한 나머지 map파일 항목 만들기 #####

            map_CHR = ['6' for i in range(0, len(map_LABELS))]
            map_GENETIC_DISTANCE = ['0' for i in range(0, len(map_LABELS))]


            df_OUTPUT_map = pd.DataFrame.from_dict({"Chr" : map_CHR, "Name" : map_LABELS.tolist(), "GD" : map_GENETIC_DISTANCE, "POS" : map_POS}).loc[:, ["Chr", "Name", "GD", "POS"]]
            print(std_MAIN_PROCESS_NAME + "Output .map file.\n")
            print(df_OUTPUT_map.head(50))

            df_OUTPUT_map.to_csv('.'.join([_OUTPUT,'HLA.map']), sep='\t', header=False, index=False)



        if MAKING_OUTPUT_PED:

            ########## < 4. Making OUTPUT.ped file > ##########

            print(std_MAIN_PROCESS_NAME + "[4] Making .ped file.\n")


            """
            
                    to_df_OUTPUT_ped = []
            
                    for i in range(0, len(INPUT_PED.index)):
            
                        line = tuple(INPUT_PED.iloc[i, :])
            
                        to_df_OUTPUT_ped.extend([
                            line[0:6] +
                            PrintGenotypes("A", line[6], line[7], sorted_keys) +
                            PrintGenotypes("C", line[10], line[11], sorted_keys) +
                            PrintGenotypes("B", line[8], line[9], sorted_keys) +
                            PrintGenotypes("DRB1", line[20], line[21], sorted_keys) +
                            PrintGenotypes("DQA1", line[16], line[17], sorted_keys) +
                            PrintGenotypes("DQB1", line[18], line[19], sorted_keys) +
                            PrintGenotypes("DPA1", line[12], line[13], sorted_keys) +
                            PrintGenotypes("DPB1", line[14], line[15], sorted_keys)
                        ])
            
                    # 실제로 HLA allele들이 genome상에서 저 position순으로 존재하기 때문에 저렇게 구성시키는게 맞다고함
                    # 가지고 시작하는 HAPMAP_CEU_HLA.ped파일의 HLA allele순서 구성은 정말 우리들끼리 임의대로 그렇게 구성해놓은 것일 뿐임.
                    # from 승호쌤's comment
                    
                    df_OUTPUT_ped = pd.DataFrame(to_df_OUTPUT_ped)
                    df_OUTPUT_ped.to_csv(_OUTPUT + '.HLA.ped', sep='\t', header=False, index=False)
                        
                    """


            to_df_OUTPUT_ped = []

            # for i in range(0, 5):
            for i in range(0, INPUT_PED.shape[0]):

                # print("\n================\n")

                line_INPUT_PED = tuple(INPUT_PED.iloc[i, :])
                # print(line_INPUT_PED)


                t_line_OUTPUT_PED = [PrintGenotypes3(line_INPUT_PED[2*j], line_INPUT_PED[2*j+1], dict_ALL_ALLELES[HLA_names[j]]) for j in range(0, len(HLA_names))]
                # print(t_line_OUTPUT_PED)
                # print(pd.Series(dict_ALL_ALLELES["A"]))

                # Flattening
                line_OUTPUT_PED = [item for eachlist in t_line_OUTPUT_PED for item in eachlist]

                # print("\nFlattened t_line_OUTPUT_PED is \n")
                # print(line_OUTPUT_PED)

                to_df_OUTPUT_ped.append(line_OUTPUT_PED)


            df_OUTPUT_ped = pd.DataFrame(to_df_OUTPUT_ped)
            df_OUTPUT_ped.index = INPUT_PED.index

            print(df_OUTPUT_ped.head())


            df_OUTPUT_ped.to_csv('.'.join([_OUTPUT, 'HLA.ped']), sep='\t', header=False, index=True)

    else:

        ### Acquiring `HLA_allele_sets_4field`.

        HLA_allele_sets_4field = {HLA_names[i]: [] for i in range(0, len(HLA_names))}

        with open(_CHPED, 'r') as f_chped:

            count = 0

            for l in f_chped:

                """
                l[:6] := ("FID", "IID", "PID", "MID", "Sex", "Phe")
                l[6:8] := HLA-A
                l[8:10] := HLA-B
                ...
                l[20:22] := HLA-DRB1
                """

                t_line = re.split(r'\s+', l.rstrip('\n'))

                for i in range(0, len(HLA_names)):

                    idx1 = 2*i + 6
                    idx2 = idx1 + 1

                    if t_line[idx1] != "0" and (t_line[idx1] not in HLA_allele_sets_4field[HLA_names[i]]):
                        HLA_allele_sets_4field[HLA_names[i]].append(t_line[idx1])
                    if t_line[idx2] != "0" and (t_line[idx2] not in HLA_allele_sets_4field[HLA_names[i]]):
                        HLA_allele_sets_4field[HLA_names[i]].append(t_line[idx2])


                count += 1
                # if count > 5 : break


        # Result checking
        print("\n4-field HLA alleles.")
        for k, v in HLA_allele_sets_4field.items():
            print("{}: {}".format(k, v))



        ### Acquiring `HLA_allele_sets_1field`

        p = re.compile(r'(\w+\*\d{2,3})') if not __previous_version else re.compile(r'(\w+:\d{2})')

        HLA_allele_sets_1field = {HLA_names[i]: [] for i in range(0, len(HLA_names))}

        for i in range(0, len(HLA_names)):

            for al_4field in HLA_allele_sets_4field[HLA_names[i]]:

                m = p.match(al_4field)

                if m and (m.group() not in HLA_allele_sets_1field[HLA_names[i]]) and (m.group() not in HLA_allele_sets_4field[HLA_names[i]]):
                    HLA_allele_sets_1field[HLA_names[i]].append(m.group())

        # Result checking
        print("\n1-field HLA alleles.")
        for k, v in HLA_allele_sets_1field.items():
            print("{}: {}".format(k, v))



        ### Merging 1field and 4field alleles

        HLA_allele_sets_merged = {}

        for i in range(0, len(HLA_names)):
            HLA_allele_sets_merged[HLA_names[i]] = HLA_allele_sets_4field[HLA_names[i]]
            HLA_allele_sets_merged[HLA_names[i]].extend(HLA_allele_sets_1field[HLA_names[i]])

            # sorting
            HLA_allele_sets_merged[HLA_names[i]].sort()

        # Result checking
        print("\nMerged and sorted.")
        for k, v in HLA_allele_sets_merged.items():
            print("{}: {}".format(k, v))



        ### Making new *.HLA.map file.

        # Concatening `HLA_allele_sets_merged` into a single list.
        l_HLA_allele_sets_merged = []
        l_HLA = []

        for i in (range(0, len(HLA_names)) if not __previous_version else HLA_names2):
            l_HLA_allele_sets_merged.extend(HLA_allele_sets_merged[HLA_names[i]])
            l_HLA.extend([HLA_names[i] for z in range(0, len(HLA_allele_sets_merged[HLA_names[i]]))])

        print("\nMerged as list.")
        print(l_HLA_allele_sets_merged)
        print(l_HLA)


        if __previous_version:
            l_HLA_allele_sets_merged = list(map(lambda x : re.sub(':', repl='', string=re.sub(':', repl='_', count=1, string=x)), l_HLA_allele_sets_merged))


        map_LABELS = ["HLA_" + al for al in l_HLA_allele_sets_merged]
        map_POS = [str(genepos_hg[_hg][_hla]) for _hla in l_HLA]

        with open(_OUTPUT + ".map", 'w') as f_HLA_map:
            f_HLA_map.writelines(('\t'.join(["6", map_LABELS[i], "0", map_POS[i]]) + "\n" for i in range(0, len(map_LABELS))))

            if __addDummyMarker:
                f_HLA_map.write('\t'.join(["6", "dummy_marker", "0", "33999999"]) + "\n")



        ### Making a new *.HLA.ped file.

        with open(_OUTPUT + ".ped", 'w') as f_HLA_ped:
            f_HLA_ped.writelines(MakeHLAPed(_CHPED, HLA_allele_sets_merged,
                                            __asCapital=__asCapital, __addDummyMarker=__addDummyMarker,
                                            __previous_version=__previous_version))




    return 0


def PrintGenotypes3(_allele1, _allele2, _seg_ALL_ALLELES):

    l_output = []

    # print("\nAlleles : {0} and {1}".format(_allele1, _allele2))
    # print("\ndict_ALL_ALLELES: \n{0}".format(_seg_ALL_ALLELES))

    if len(_seg_ALL_ALLELES) > 0:

        if ((_allele1 != "0") and (_allele2 != "0")):  # The condition for checking integer value 0 won't be included here because .ped file was read with "dtype=str" option.

            t_sr1 = pd.Series(_seg_ALL_ALLELES, index=pd.Index(_seg_ALL_ALLELES)).apply(lambda x: "P" if (x in _allele1) else "A")
            # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr1))
            t_sr2 = pd.Series(_seg_ALL_ALLELES, index=pd.Index(_seg_ALL_ALLELES)).apply(lambda x: "P" if (x in _allele2) else "A")
            # print("{0} : t_sr1 is \n{1}".format(_locus, t_sr2))

            for i in range(0, len(t_sr1)):
                l_output.append(t_sr1.iat[i])
                l_output.append(t_sr2.iat[i])

            # print(l_output)
            return l_output


        else:

            # Comments... at most below part.
            # print("At least one of allele is 0")

            for i in range(0, len(_seg_ALL_ALLELES)):
                l_output.append("0")
                l_output.append("0")

            return l_output

    else:
        # In cases such as "DPA1" or "DPB1 where any alleles don't appear, Just skip.
        # Then just return NULL
        return l_output


# (2019. 1. 3.) Introduced for memory issues.
def PrintGenotypes4(_allele1, _allele2, _HLA_allele_sets_merged_byHLA, __asCapital=True):

    l_output = []

    _present_ = "P" if __asCapital else "p"
    _absent_ = "A" if __asCapital else "a"

    if len(_HLA_allele_sets_merged_byHLA) > 0:

        if ((_allele1 != "0") and (_allele2 != "0")):

            allele1_PAcheck = [_present_ if (x in _allele1) else _absent_ for x in _HLA_allele_sets_merged_byHLA]
            allele2_PAcheck = [_present_ if (x in _allele2) else _absent_ for x in _HLA_allele_sets_merged_byHLA]

            # "in" operator를 쓰면 안될거같은게, 예를 들어 "10"은 0"10"1 이거 이런식으로 찾아와서 P박아버리는 케이스가 좀 있을거같음.

            for i in range(0, len(allele1_PAcheck)):
                l_output.append(allele1_PAcheck[i])
                l_output.append(allele2_PAcheck[i])

        else:

            for i in range(0, len(_HLA_allele_sets_merged_byHLA)):
                l_output.append("0")
                l_output.append("0")


        return '\t'.join(l_output)


    else:

        # In cases such as "DPA1" or "DPB1 where any alleles don't appear, Just skip.
        pass



def MakeHLAPed(_CHPED, _HLA_allele_sets_merged, __asCapital=True, __addDummyMarker=False, __previous_version=False):

    with open(_CHPED, 'r') as f_chped:

        count = 0

        for l in f_chped:

            t_line = re.split(r'\s+', l.rstrip('\n'))

            """
            t_line[:6] := ("FID", "IID", "PID", "MID", "Sex", "Phe")
            t_line[6:8] := HLA-A
            t_line[8:10] := HLA-B
            ...
            t_line[20:22] := HLA-DRB1
            """

            __ped_info__ = '\t'.join(t_line[:6])
            __genomic_info__ = '\t'.join([
                PrintGenotypes4(t_line[2 * i + 6], t_line[2 * i + 7], _HLA_allele_sets_merged[HLA_names[i]], __asCapital=__asCapital)
                for i in (range(0, len(HLA_names)) if not __previous_version else HLA_names2) if len(_HLA_allele_sets_merged[HLA_names[i]]) > 0
            ])

            __return__ = '\t'.join([__ped_info__, __genomic_info__])


            if __addDummyMarker:
                dummy_markers = '\t'.join(['d', 'D'] if bool(count % 2) else ['D', 'd'])
                __return__ = '\t'.join([__return__, dummy_markers])


            yield __return__ + "\n"

            count += 1





if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        encodeHLA.py

        This script encodes HLA alleles into bi-allelic markers (for imputation and PLINK analysis)
        The input ped file should contain: FID,IID,pID,mID,SEX,PHENO,
                                            2 each of: HLA-A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1

    ###########################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-ped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)

    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"], metavar="hg", default="19")

    parser.add_argument("--previous-version", help="\nIf you give this option, The MakeReference will work as original version.\n\n",
                        action='store_true')



    ##### <for Test> #####

    # (2018. 2. 28.)
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./TEST_0228", "--hg", "19"])

    # args = parser.parse_args(["./HAPMAP_CEU_HLA.ped", "./TEST_0228", "--hg", "19"])
    # args = parser.parse_args(["./COATING_TEST.coated.txt", "./TEST_0305", "-hg", "19"])

    # (2018. 7. 16.)

    # args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-hg", "19",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_RECODE_v2/makereference_recode_v2/MODULE_TEST_HAPMAP_CEU_HLA.4field.imgt370"])

    # (2019. 01. 06.)

    # --previous-version
    args = parser.parse_args(["-ped", "/Users/wansun/Git_Projects/MakeReference_v2/data/MakeReference/HAPMAP_CEU_HLA.old.chped",
                              "-hg", "18",
                              "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190106/_3_encodeHLA/_Case_HAPMAP_CEU.HLA",
                              "--previous-version"])



    ##### <for Publication> #####

    # args = parser.parse_args()
    print(args)


    encodeHLA(args.ped, args.o, args.hg, __previous_version=args.previous_version)
