import os
import sys

def find_slha_directories(dirlist, filename_template=".slha"):
    counter=0
    slha_file_list = []
    for iroot in dirlist:
        print("Searching for data in {}".format(iroot))
        for root, dirs, files in os.walk(iroot, followlinks=True):
            for ifile in files:
                if (filename_template in ifile):
                    slha_file_list[len(slha_file_list):] = ["{}/{}".format(root, ifile)]
    if len(slha_file_list)>0:
        print("Found a total of {} data dir".format(len(slha_file_list)))
    else:
        print("No data found in {}".format(dirlist))
        sys.exit(1)
    return slha_file_list


def load_data(data_file_list, data_description, decay_list = []):
    data_dict = {}
    decay_found_flag = {}
    for idata in data_description:
        data_dict[idata[3]] = []
    for idecay in decay_list:
        data_dict[idecay[-1]] = []
    # Loop aver the files
    for file_counter, ifile in enumerate(data_file_list):
        if ((file_counter % 1000) == 0):
            print("At file {}".format(file_counter))
        with open(ifile, 'r') as fh:
            filetext = fh.read()
        current_block = None
        decay_block_flg = False
        for idecay in decay_list:
            decay_found_flag[idecay[-1]] = False
        for iline in filetext.splitlines():
            isplit = iline.split()
            if isplit[0] == '#':
                continue
            if ((isplit[0].lower() == "block") or (isplit[0].upper() == "DECAY")):
                current_block = isplit[1]
                current_block_data = isplit[2:]
            slha_comment = isplit[-1]
            value_to_extract = list(filter(lambda x: (x[0] == current_block and x[1] == slha_comment) ,
                                                      data_description))
            if len(value_to_extract) > 0:
                if len(value_to_extract) > 1:
                    print("Error: more than matching for line {}".format(isplit))
                    print("Matching data definitions are {}".format(value_to_extract))
                data_key = value_to_extract[0][3]
                data_slha_idx = value_to_extract[0][2]
                data_dict[data_key].append(float(isplit[data_slha_idx]))
            if isplit[0].upper() == "DECAY":
                decay_block_flg = False
                # Is it a decay block we are searching for?
#                print(isplit)
                # Match if the particle decaying is specific in the DECAY entry
                decay_to_extract_list = list(filter(lambda x: (x[0] == float(isplit[1])), decay_list))
                if len(decay_to_extract_list) > 0:
                    decay_block_flg = True
            if (decay_block_flg == True):
                try:
                    # 2->2 decays, check the two decay products
                    if isplit[1] == "2":
                        matched_decay = list(filter(lambda x: (x[1] == float(isplit[2]) and x[2] == float(isplit[3])),
                                                    decay_to_extract_list))
                    # 2->3 decays, check the two decay products
                    elif isplit[1] == "3":
                        matched_decay = list(filter(lambda x: (x[1] == float(isplit[2]) and x[2] == float(isplit[3]) and x[2] == float(isplit[3])),
                                                    decay_to_extract_list))
                    else:
                        matched_decay = []
                except ValueError:
                    pass
                else:
                    if (len(matched_decay)>0):
 #                       print("matched decay", matched_decay)
 #                       print(isplit)
                        # The BR is always the first entry; in the decay definition the last element of the list is the dictionary key to be used
                        data_dict[matched_decay[0][-1]].append(float(isplit[0]))
                        decay_found_flag[matched_decay[0][-1]] = True
        for idecay in decay_list:
            if decay_found_flag[idecay[-1]] == False:
#                print("BR {} set to zero".format(idecay))
                data_dict[idecay[-1]].append(0.)
#        sys.exit(0)
    print("Sanity check")
    nfiles = len(data_file_list)
    print("Read a total of {} files".format(nfiles))
    for ikey in data_dict:
        print("Len data_dict[{}] = {}".format(ikey, len(data_dict[ikey])))
        assert(len(data_dict[ikey]) == nfiles)
    return data_dict
