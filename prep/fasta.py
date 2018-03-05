"""
FASTA File Processing
"""

import json
import glob
import re

from misc.string import *


def flybase_fasta_to_dict(infile):
    """
    Convert FASTA to Python dictionary
    :param infile: FASTA file from FlyBase
    :return: Python dictionary
    """

    fasta_dict = dict()
    header = None

    fin = open(infile, 'r')
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

            fasta_dict[header] = entry_dict

        else:

            if "sequence" not in fasta_dict[header].keys():
                fasta_dict[header]["sequence"] = []
            fasta_dict[header]["sequence"].append(line.strip())

    fin.close()

    return fasta_dict


def wormbase_fasta_to_dict(infile):
    """
    Convert FASTA to Python dictionary
    :param infile: FASTA file from WormBase
    :return: Python dictionary
    """

    fasta_dict = dict()
    header = None

    fin = open(infile, 'r')
    for line in fin.readlines():

        if line.startswith(">"):

            entry_dict = dict()

            header = re.findall(r">(\S+)", line)[0]

            search_entry = re.findall(r"(\S+)=([^\s\"]+)", line) + re.findall(r"(\S+)=\"(.+?)\"", line)
            for item in search_entry:
                key, value = item
                entry_dict[key] = value

            fasta_dict[header] = entry_dict

        else:

            if "sequence" not in fasta_dict[header].keys():
                fasta_dict[header]["sequence"] = []
            fasta_dict[header]["sequence"].append(line.strip())

    fin.close()

    return fasta_dict


def fasta_to_json(infile, outfile, db_name):
    """
    Convert FASTA to JSON
    :param db_name: DataBase source, e.g., FlyBase, WormBase
    :param infile: FASTA file
    :param outfile: JSON file
    :return:
    """

    if db_name == "flybase":
        json_dict = flybase_fasta_to_dict(infile)
    elif db_name == "wormbase":
        json_dict = wormbase_fasta_to_dict(infile)

    fout = open(outfile, 'w')
    json.dump(json_dict, fout, indent=5)
    fout.close()


def fasta_to_json_folder(infolder, outfolder, db_name):
    """
    Convert FASTA to JSON
    :param db_name: DataBase source, e.g., FlyBase, WormBase
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

        fasta_to_json(file, outfile, db_name)


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

            header = re.findall(r"(>\S+)", line)[0]
            print(header, file=fout)

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


def fasta_to_conversion_dict(infile, db_name):
    """
    Convert FASTA file to Python dictionary containing Protein ID -> Gene ID
    :param infile: FASTA file
    :param db_name: DataBase source, e.g., FlyBase, WormBase
    :return:
    """

    id_conversion_dict = dict()

    if db_name == "wormbase":
        fasta_dict = wormbase_fasta_to_dict(infile)
        for key, value in fasta_dict.items():
            id_conversion_dict[key] = value["gene"]
    elif db_name == "flybase":
        fasta_dict = flybase_fasta_to_dict(infile)
        for key, value in fasta_dict.items():
            id_conversion_dict[key] = value["parent"][0]

    return id_conversion_dict


def fasta_to_conversion_file(infile, outfile, db_name):
    """
    Convert FASTA file to conversion file containing Protein ID -> Gene ID
    :param db_name: DataBase source, e.g., FlyBase, WormBase
    :param outfile: ID Conversion file in .tsv format
    :param infile: FASTA file
    :return:
    """

    fout = open(outfile, "w")

    id_conversion_dict = fasta_to_conversion_dict(infile, db_name)
    for key, value in id_conversion_dict.items():
        print(key, value, sep="\t", file=fout)

    fout.close()


def fasta_to_conversion_folder(infolder, outfolder, db_name):
    """
    Convert FASTA file to ID Conversion file containing Protein ID -> Gene ID
    :param db_name: DataBase source, e.g., FlyBase, WormBase
    :param infolder: FASTA folder
    :param outfolder: ID Conversion folder
    :return:
    """

    infolder = check_folder_path(infolder)
    outfolder = check_folder_path(outfolder, True)

    for file in glob.glob(infolder + '*'):
        filename = os.path.basename(file)
        new_filename = filename.split(".")[0] + ".tsv"
        outfile = outfolder + new_filename

        fasta_to_conversion_file(file, outfile, db_name)

