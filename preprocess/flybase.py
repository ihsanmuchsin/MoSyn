"""
Pre processing data from FlyBase
"""

import json
import glob
import re

from misc.string import *


def fasta_to_json(infile, outfile):
    """
    Convert FASTA to JSON
    :param infile: FASTA file
    :param outfile: JSON file
    :return:
    """

    fin = open(infile, 'r')

    json_dict = dict()
    header = None

    for line in fin.readlines():

        if line.startswith(">"):

            entry_dict = dict()

            header = re.findall(r">(\S+)", line)[0]

            search_entry = re.findall(r"(\S+)=(\S+);", line)
            for item in search_entry:
                key, value = item

                if key == 'loc':
                    subentry = re.findall(r"([^\s,]+):(\S+\)+)", value)
                else:
                    subentry = re.findall(r"([^\s,]+):([^\s,]+)", value)

                if subentry:
                    value = dict()
                    for itm in subentry:
                        k, v = itm
                        value[k] = v

                if type(value) == str and len(value.split(",")) > 1:
                    value = value.split(",")

                entry_dict[key] = value

            json_dict[header] = entry_dict

        else:

            if "sequence" not in json_dict[header].keys():
                json_dict[header]["sequence"] = []
            json_dict[header]["sequence"].append(line.strip())

    fin.close()

    fout = open(outfile, 'w')
    json.dump(json_dict, fout, indent=5)
    fout.close()


def fasta_to_json_folder(infolder, outfolder):
    """
    Convert FASTA to JSON
    :param infolder: FASTA folder
    :param outfolder: JSON folder
    :return:
    """

    infolder = check_folder_path(infolder)
    outfolder = check_folder_path(outfolder, True)

    for file in glob.glob(infolder + '*'):
        filename = os.path.basename(file)
        new_filename = filename.split(".")[0] + ".json"
        outfile = outfolder + new_filename

        fasta_to_json(file, outfile)


def clean_header_fasta(infile, outfile):
    """
    Clean FASTA header
    :param infile: FASTA file
    :param outfile: FASTA file with clean header
    :return:
    """

    fout = open(outfile, 'w')
    fin = open(infile, 'r')

    for line in fin.readlines():
        if line.startswith(">"):

            line_elem = line.strip().split(" ")
            print(line_elem[0], file=fout)

        else:
            fout.write(line)

    fin.close()
    fout.close()


def clean_header_fasta_folder(infolder, outfolder):
    """
    Clean FASTA header
    :param infolder: FASTA folder
    :param outfolder: FASTA with clean header folder
    :return:
    """

    infolder = check_folder_path(infolder)
    outfolder = check_folder_path(outfolder, True)

    for file in glob.glob(infolder + '*'):
        filename = os.path.basename(file)
        outfile = outfolder + filename

        clean_header_fasta(file, outfile)
