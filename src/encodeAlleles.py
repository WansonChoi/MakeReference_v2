import os, sys, re
import argparse, textwrap


# Global Variables
Normal_Bases = ['A', 'T', 'C', 'G']

def encodeAllele(*args, **kwargs):


    _tped, _out = args
    _function_switch = "encoding" if bool(kwargs["_encode"]) else "decoding"


    p = re.compile(r'\.tped$')
    _out = _out if not p.search(_out) else re.sub(p, '', _out)



    if _function_switch == "encoding":

        # print("\n[Encoding]")

        __EncodeTable__ = getEncodeTable(_tped)

        # for k,v in __EncodeTable__.items():
        #     print("{} : {}".format(k,v))


        ### Main encoding
        # print("Processing ENcoding *.tped file.")
        with open(_out + ".aENCODED.tped", 'w') as f_out:
            f_out.writelines(EncodeAllele(_tped, __EncodeTable__, _sep=' '))


        ### Generate Decoding Table.
        # print("Generating *.emap file.")
        with open(_out + ".aENCODED.emap", 'w') as f_dmap:
            f_dmap.writelines((' : '.join([k, str(v)])+"\n" for k, v in __EncodeTable__.items()))


        # print("\n<Result(s)>")
        # print("(1) Encoded *.tped file : {}".format(_out + ".aENCODED.tped"))
        # print("(2) Encoding rule table : {}".format(_out + ".aENCODED.emap"))


        return [_out + ".aENCODED.tped", _out + ".aENCODED.emap"]



    elif _function_switch == "decoding":

        # print("\n[Decoding]")

        if not bool(kwargs["_emap"]):
            print("Error. The argument \"-emap\" must be given.\n")
            sys.exit()
        else:
            _emap = kwargs["_emap"]

        # print("_emap : {}".format(_emap))


        ### Retrieving Factors
        # print("Getting DEcoding table.")
        __DecodeTable__ = getDecodeTable(_emap)

        # for k, v in __DecodeTable__.items():
        #     print("{} : {}".format(k, v))


        ### Main Decoding
        # print("Processing decoding *.tped file.")
        with open(_out + ".aDECODED.tped", 'w') as f_out:
            f_out.writelines(DecodeAllele(_tped, __DecodeTable__, _sep=' '))



        # print("\n<Result(s)>")
        # print("Decoded *.tped file : {}".format(_out + ".aDECODED.tped"))



        return _out + ".aDECODED.tped"





def getEncodeTable(_tped):

    __EncodeTable__ = {}
    count = 0

    with open(_tped, 'r') as f_tped:

        for l in f_tped:

            line = re.split(r'\s+', l.rstrip('\n'))

            line_MarkerInfo = line[:4]
            line_GenomicInfo = line[4:]

            t_line_GenomicInfo = line_GenomicInfo
            t_line_GenomicInfo.sort()
            t_factors = set(t_line_GenomicInfo).difference({"0"}) # Exclude "0"
            # print("{} : {}".format(line_MarkerInfo, t_factors))


            if len(t_factors) == 2:

                if (t_factors == {"A", "C"} or t_factors == {"A", "G"} or t_factors == {"A", "T"} or
                        t_factors == {"C", "G"} or t_factors == {"C", "T"} or t_factors == {"G", "T"}):

                    # 6 cases(choose(4,2)) of normal allele combinations.
                    pass

                else:
                    # The cases out of this must be re-encoded
                    # print("The factor set {} will be re-encoded.".format(t_factors))

                    t_factors = list(t_factors)
                    t_factors.sort()

                    __EncodeTable__[line_MarkerInfo[1]] = {t_factors[i]: Normal_Bases[i] for i in range(0, len(t_factors))}

                    # print("It is re-encoded into {}".format(__EncodeTable__[line_MarkerInfo[1]]))


            count += 1
            # if count > 10 : break;


    return __EncodeTable__



def getDecodeTable(_emap):

    __DecodeTable__ = {}

    p = re.compile(r'\s+:\s+')
    count = 0

    with open(_emap, 'r') as f_emap:

        for l in f_emap:

            line = re.split(p, l.rstrip('\n'))
            # print(line)


            t_dict = eval(line[1])
            __DecodeTable__[line[0]] = {v: k for k, v in t_dict.items()}


            # print("{} : {}\n".format(line[0], __DecodeTable__[line[0]]))


            count += 1
            # if count > 50 : break



    return __DecodeTable__




def EncodeAllele(_tped, _EncodeTable, _sep='\t'):

    Markers = _EncodeTable.keys()
    # print(Markers)
    count = 0

    with open(_tped, 'r') as f_tped:

        for l in f_tped:

            line = re.split(r'\s+', l.rstrip('\n'))
            line_MarkerInfo = line[:4]


            if line_MarkerInfo[1] in Markers:

                t_keys = _EncodeTable[line_MarkerInfo[1]].keys()

                line_GenomicInfo = line[4:]

                for i in range(0, len(line_GenomicInfo)):

                    if line_GenomicInfo[i] in t_keys:
                        line_GenomicInfo[i] = _EncodeTable[line_MarkerInfo[1]][line_GenomicInfo[i]]


                yield _sep.join([_sep.join(line_MarkerInfo), _sep.join(line_GenomicInfo)]) + "\n"

            else:
                yield l



            count += 1
            # if count > 10 : break




def DecodeAllele(_tped, _DecodeTable, _sep='\t'):

    Markers = _DecodeTable.keys()
    # print(Markers)
    count = 0

    with open(_tped, 'r') as f_tped:
        for l in f_tped:

            line = re.split(r'\s+', l.rstrip('\n'))
            line_MarkerInfo = line[:4]


            if line_MarkerInfo[1] in Markers:

                t_keys = _DecodeTable[line_MarkerInfo[1]].keys()

                line_GenomicInfo = line[4:]

                for i in range(0, len(line_GenomicInfo)):

                    if line_GenomicInfo[i] in t_keys:
                        line_GenomicInfo[i] = _DecodeTable[line_MarkerInfo[1]][line_GenomicInfo[i]]


                yield _sep.join([_sep.join(line_MarkerInfo), _sep.join(line_GenomicInfo)]) + "\n"

            else:
                yield l



            count += 1
            # if count > 5 : break







if __name__ == "__main__" :

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        encodeAllele.py        
        
        - Glue module to handle compatibility issues related to beagle4.x
        - In beagle4.x(> v3.x.x), 
            (1) No other allele characters execpt 'A', 'C', 'G', 'T' and 'N' are allowed,
            (2) No same genomic position is allowed to each markers.
        - In beagle v3.x.x, above two issues weren't dealt with at all, which are the reasons why
            S. Jia and B. Han chose beagle(v3.x.x) as main engine for 'MakeReference' and 'SNP2HLA'.
        - By the way, This module solves the first one in the above 2 issues by encoding {'P', 'A'},
            {'S', 'L'} or other allele combinations into the {'A', 'T'}.
        - The table containing encoding rules are generated with encoding job 
            so that decoding them back to original values can be done.


    #################################################################################################
                                     '''),
                                     add_help=False)


    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    func_switch = parser.add_mutually_exclusive_group(required=True)
    func_switch.add_argument("--encode", help="\nEncode *.tped file.\n\n", action='store_true')
    func_switch.add_argument("--decode", help="\nDecode *.tped file.\n\n", action='store_true')

    parser.add_argument("-tped", help="\nTransposed genotype data(*.tped)\n\n", required=True)
    parser.add_argument("-emap", help="\nThe one of outputs(*.emap) after encoding a *.tped file.\n\n")
    parser.add_argument("-o", help="\nOutput file prefix\n\n", required=True)

    # parser.add_argument("--allele-before", help="\nTarget allele character in the given *.ped file.\n\n", required=True)
    # parser.add_argument("--allele-after", help="\nThe new allele character which an user wants to change.\n\n", required=True)





    ##### <for Test> #####

    # Partial Test Data
    # args = parser.parse_args(["--encode",
    #                           "-tped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190129/20190129_Panel.tped",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190129/20190129_Panel.bglv4"
    #                           ])

    # Test Data
    # args = parser.parse_args(["--decode",
    #                           "-tped", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190129/20190129_Panel.bglv4.tped",
    #                           "-emap", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190129/20190129_Panel.bglv4.emap",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190129/20190129_Panel.bglv4.decoded"
    #                           ])


    ##### <for Publication> #####

    args = parser.parse_args()
    print(args)




    # Summary print
    print("\n<Argument Checking>")
    print("_tped : {}".format(args.tped))
    print("_out : {}".format(args._out))
    _function_switch = "encoding" if args.encode else "decoding"
    print("_switch : {}".format(_function_switch))



    encodeAllele(args.tped, args.o, _encode=args.encode, _decode=args.decode, _emap=args.emap)

