"""
GTF File Processing
"""

import json
import glob
import re
import pandas as pd
import os

from misc.string import check_folder_path


def gtf_to_dict(infile):
    """
    Convert GTF to python dictionary
    :param infile: GTF file
    :return: Python dictionary
    """

    gtf_dict = dict()

    fin = open(infile, 'r')
    for line in fin.readlines():

        if not line.startswith("#"):

            entry_dict = dict()

            line_elem = line.strip().split("\t")
            seqname = line_elem[0]

            entry_dict["source"], entry_dict["feature"], entry_dict["start"], entry_dict["end"], entry_dict["score"], \
            entry_dict["strand"], entry_dict["frame"] = line_elem[1:-1]

            attribute = dict()
            attribute_items = re.findall(r"(\S+)\s\"(\S+)\";", line_elem[-1])
            for item in attribute_items:
                key, value = item
                attribute[key] = value
            entry_dict["attribute"] = attribute

            if seqname not in gtf_dict.keys():
                gtf_dict[seqname] = []
            gtf_dict[seqname].append(entry_dict)

    fin.close()

    return gtf_dict


def gtf_to_json(infile, outfile):
    """
    Convert GTF to JSON
    :param infile: GTF file
    :param outfile: JSON file
    :return:
    """
    json_dict = gtf_to_dict(infile)

    fout = open(outfile, 'w')
    json.dump(json_dict, fout, indent=5)
    fout.close()


def gtf_to_json_folder(infolder, outfolder):
    """
    Convert GTF to JSON
    :param infolder: GTF folder
    :param outfolder: JSON folder
    :return:
    """

    infolder = check_folder_path(infolder)
    outfolder = check_folder_path(outfolder, True)

    for file in glob.glob(infolder + '*'):
        filename = os.path.basename(file)
        new_filename = filename.split(".")[0] + ".json"
        outfile = outfolder + new_filename

        gtf_to_json(file, outfile)


def gtf_to_iadhore_list(infile, outfolder):
    """
    Convert GTF file to i-ADHoRe list
    :param infile: GTF file
    :param outfolder: i-ADHoRe list
    :return:
    """

    outfolder = check_folder_path(outfolder, True)
    gtf_dict = gtf_to_dict(infile)
    suboutfolder = outfolder + os.path.basename(infile).split('.')[0]
    suboutfolder = check_folder_path(suboutfolder, True)

    # Create pandas dataframe
    gene_location_list = []
    for key, values in gtf_dict.items():
        for value in values:
            if value['feature'] == 'gene':
                gene_location_list.append((value['attribute']['gene_id'], value['strand'],
                                           key, int(value['start']), int(value['end'])))
    labels = ['gene', 'strand', 'chromosome', 'start', 'end']
    df = pd.DataFrame.from_records(gene_location_list, columns=labels)
    groups = df.groupby(['chromosome'])

    for chromosome, group in groups:
        outfile = suboutfolder + chromosome + '.lst'
        group_sort = group.sort_values(['start'])

        fout = open(outfile, 'w')
        for idx in group_sort.index:
            print(group_sort.at[idx, 'gene'] + group_sort.at[idx, 'strand'], file=fout)
        fout.close()


def gtf_to_iadhore_list_folder(infolder, outfolder):
    """
    Convert GTF file to i-ADHoRe list
    :param infolder: GTF folder
    :param outfolder: i-ADHoRe list
    :return:
    """

    infolder = check_folder_path(infolder)
    outfolder = check_folder_path(outfolder, True)

    for file in glob.glob(infolder + '*'):
        gtf_to_iadhore_list(file, outfolder)


def gtf_to_dict_folder(infolder):
    """
    Convert GTF files to dictionary
    :param infolder: GTF folder
    :return: Python dictionary
    """

    all_gtf = dict()

    infolder = check_folder_path(infolder)
    for file in glob.glob(infolder + '*'):
        filename = os.path.basename(file)
        filename = filename.split(".")[0]

        all_gtf[filename] = gtf_to_dict(file)

    return all_gtf

