"""
Pre processing data from WormBase
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
    key_temp = 0

    for line in fin.readlines():

        if line.startswith(">"):

            test0 = re.findall(r">(\S+)", line)
            print(test0)
            test1 = re.findall(r"(\w+)=\"(.+?)\"", line)
            print(test1)
            test2 = re.findall(r"(\S+)=([^\s_\"]+)", line)
            print(test2)

            key_temp += 1

            line_elem = line.strip().split(" ")

            entry_dict = dict()
            protein_id = line_elem[0].strip(">")
            entry_dict["protein_id"] = protein_id
            print(protein_id)

            if len(line_elem) > 1:
                for i in range(1, len(line_elem)):
                    print(line_elem[i])
                    key, value = line_elem[i].split('=')
                    entry_dict[key] = value

            json_dict[key_temp] = entry_dict

        else:

            if "sequence" not in json_dict[key_temp].keys():
                json_dict[key_temp]["sequence"] = []
            json_dict[key_temp]["sequence"].append(line.strip())

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
