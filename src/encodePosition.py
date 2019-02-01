import os, sys, re
import argparse, textwrap


def encodePosition(*args, **kwargs):

    p = re.compile(r'\.bim$|\.map$')
    _out = args[0] if not p.search(args[0]) else re.sub(p, '', args[0])

    _input = kwargs["_map"] if bool(kwargs["_map"]) else kwargs["_bim"]
    f_isMap = bool(kwargs["_map"])
    _function_switch = "encoding" if bool(kwargs["_encode"]) else "decoding"

    # print("\n<Argument Checking>")
    #
    # print("_input : {}".format(_input))
    # print("_out : {}".format(_out))
    # print("_switch : {}".format(_function_switch))


    if _function_switch == "encoding":
        # print("\n[Encoding]")

        # Making a genomic position tree
        __GenPosTree__, __MarkerToPos__ = getGenDictionary(_input)


        # for k, v in __GenPosTree__.items():
        #     print("{} : {}".format(k, v))


        """
        (1) Marker to Encoded genomic position mapping
        (2) Table for decoding (*.pmap file) => Marker to (encoded position to original position)

        """

        ### Main Encoding

        # Initializaiton (a Marker to `None`)


        __isOccupied__ = list(__MarkerToPos__.values())

        count = 0

        for k, v in __GenPosTree__.items():

            t_position = k  # Genomic position
            t_markers = v # a list of markers which have same genomic position.


            if len(t_markers) > 1:

                # print(t_markers)
                # print(t_position)
                #
                # print("Offset marker : {}".format(t_markers[0]))
                # print("Offset : {}".format(t_position))

                offset = t_position + 1     # Excluding the 1st element

                for j in range(1, len(t_markers)):

                    while offset in __isOccupied__:
                        offset += 1

                    __MarkerToPos__[t_markers[j]] = offset
                    __isOccupied__.append(offset)


            count += 1
            # if count > 5 : break


        ### Writing outputs

        # (1) Genomic position encoded file (*.bglv4.map)
        with open(_out+".pENCODED.{}".format("map" if f_isMap else "bim"), 'w') as f_out:
            f_out.writelines(WriteEncodedPos(_input, __MarkerToPos__))



        # ### Results
        # # print("\n<Results>")
        # # print("(1) Genomic position encoded map file : {}".format(_out+".pENCODED.{}".format("map" if f_isMap else "bim")))
        # # print("(2) Position map table file(*.pmap) : {}".format(_out+".pENCODED.pmap"))


        return _out+".pENCODED.{}".format("map" if f_isMap else "bim")




    elif _function_switch == "decoding":
        print("\n[Decoding]")

        if not kwargs["_pmap"]:
            print("Error. \"*.pmap\" file is needed to decode *.pENCODED.{bim,map} file.")
            sys.exit()

        else:
            _pmap = kwargs["_pmap"]












def getGenDictionary(_input):

    __MarkerToPos__ = {}
    count = 0

    # Initializing `__PosToMarkers__`
    with open(_input, 'r') as f_input:
        for l in f_input:

            line = re.split(r'\s+', l.rstrip('\n'))

            __MarkerToPos__[line[1]] = int(line[3])

            count += 1
            # if count > 5 : break

    __PosToMarkers__ = {pos: [] for pos in __MarkerToPos__.values()}    # Initializing `__PosToMarkers__` as null lists.


    # Clustering Marker Labels based on genomic positions.

    count = 0
    with open(_input, 'r') as f_input:
        for l in f_input:

            line = re.split(r'\s+', l.rstrip('\n'))

            t_marker = line[1]
            t_genomic_position = int(line[3])

            # print("Marker : {}".format(t_marker))
            # print("GenPos : {}".format(t_genomic_position))
            # print("list : {}\n".format(__PosToMarkers__[t_genomic_position]))

            __PosToMarkers__[t_genomic_position].append(t_marker)


            count += 1
            # if count > 5 : break


    return __PosToMarkers__, __MarkerToPos__



def WriteEncodedPos(_input, _MarkerToPos_, _sep='\t'):

    with open(_input, 'r') as f_input:
        for l in f_input:

            line = re.split(r'\s+', l.rstrip('\n'))

            t_marker = line[1]

            yield _sep.join([line[0], t_marker, line[2], str(_MarkerToPos_[t_marker]), _sep.join(line[4:])]) + "\n"


# def WritePMAP(_)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        encodePosition.py        

        - Glue module to handle compatibility issues related to beagle4.x
        - In beagle4.x(> v3.x.x), 
            (1) No other allele characters execpt 'A', 'C', 'G', 'T' and 'N' are allowed,
            (2) No same genomic position is allowed to each markers.
        - In beagle v3.x.x, above two issues weren't dealt with at all, which are the reasons why
            S. Jia and B. Han chose beagle(v3.x.x) as main engine for MakeReference and SNP2HLA.
        - By the way, This module solves the second one in the above 2 issues by encoding same genomic 
            positions to unique integer values.
        - The table containing encoding rules are generated with encoding job 
            so that decoding them back to the original values can be done.


    #################################################################################################
                                     '''),
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    func_switch = parser.add_mutually_exclusive_group(required=True)
    func_switch.add_argument("--encode", help="\nEncode *.map(or bim) file.\n\n", action='store_true')
    func_switch.add_argument("--decode", help="\nDecode *.map(or bim) file.\n\n", action='store_true')

    MaporBim = parser.add_mutually_exclusive_group(required=True)
    MaporBim.add_argument("-map", help="\nPlink *.map file.\n\n")
    MaporBim.add_argument("-bim", help="\nPlink *.bim file.\n\n")

    parser.add_argument("-pmap", help="\nEncoding rule information which is generated in encoding a *.map(or bim) file.\n\n")
    parser.add_argument("-o", help="\nOutput file prefix\n\n", required=True)


    ##### <for Test> #####

    # Partial Test Data
    # args = parser.parse_args(["--encode",
    #                           "-bim", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190130/20190129_Panel.bim",
    #                           "-o", "/Users/wansun/Git_Projects/MakeReference_v2/tests/20190130/20190201_Panel.bim"
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

    encodePosition(args.o, _encode=args.encode, _decode=args.decode, _map=args.map, _bim=args.bim, _pmap=args.pmap)

