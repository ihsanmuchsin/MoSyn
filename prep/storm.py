"""
STORM processing
"""

import glob
import os
import re

from misc.string import check_folder_path


def transfac_to_gtf(infile, outfile, binding_site_id="BS", id_start_index=0, motif_name="MOTIF", return_index=False):
    """
    Convert a TRANSFAC file / STORM output to GTF
    :param motif_name: The motif name, e.g., CTCF.
    :param binding_site_id: The ID of the binding sites, e.g., BS. It will be BS0, BS1, .., BSx in the output
    :param id_start_index: the start index of the motif ID, e.g., 0,1,2,...,x
    :param infile: TRANSFAC file / STORM output
    :param outfile: GTF file of binding sites (BS) location
    :return:
    """

    fout = open(outfile, 'w')

    fin = open(infile, 'r')
    for line in fin.readlines():

        if line.startswith("BS"):

            info = line.strip().split(";")

            sequence = info[0].split(" ")[-1].strip()
            start = int(info[2].strip())
            length = int(info[3].strip())
            gap = 0 if " " else int(info[4].strip())
            end = start + length + gap
            score = float(info[-1].strip())
            strand = "+" if info[-2].strip() == "p" else "-"

            attr = ["motif_id", "\""+binding_site_id+str(id_start_index)+"\";",
                    "motif_name", "\""+motif_name+"\";",
                    "sequence", "\""+sequence+"\";"]
            attr = " ".join(attr)

            print(info[1], "STORM", "Motif", start, end, score, strand, ".", attr, sep="\t", file=fout)

            id_start_index += 1

    fout.close()

    if return_index:
        return id_start_index


def transfac_to_gtf_folder(infolder, outfolder, binding_site_id="BS", id_start_index=0, motif_name="MOTIF"):
    """
    Convert a TRANSFAC folder / STORM output to GTF
    :param infolder: TRANSFAC / STORM Output folder
    :param outfolder: GTF folder
    :param binding_site_id: The ID of the binding sites, e.g., BS. It will be BS0, BS1, .., BSx in the output
    :param id_start_index: the start index of the motif ID, e.g., 0,1,2,...,x
    :param motif_name: The motif name, e.g., CTCF.
    :return:
    """

    infolder = check_folder_path(infolder)
    outfolder = check_folder_path(outfolder, True)

    for file in glob.glob(infolder + '*'):
        filename = os.path.basename(file)
        new_filename = filename.split(".")[0] + ".gtf"
        outfile = outfolder + new_filename

        last_index = transfac_to_gtf(file, outfile, binding_site_id, id_start_index, motif_name, True)
        id_start_index = last_index
