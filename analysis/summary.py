"""
Generating summary
"""

import glob
import os

import prep.iadhore as pi
import core.mosyn as cm

from misc.string import check_folder_path


def generate_orthofinder_summary(infolder, outfile, genus_index=-3, alignment_index=-2):
    """
    Generate orthofinder summary.
    Structure of the folder Genus -> Alignment Method -> Result
    :param outfile: Output file
    :param alignment_index: Alignment folder index from result
    :param genus_index: Genus folder index from result
    :param infolder: Folder containing orthofinder result
    :return:
    """

    fout = open(outfile, 'w')

    print("Genus", "Alignment", "Number_of_Genes", "Number_of_Genes_in_Orthogroups", "Number_of_Unassigned genes",
          "Number_of_Orthogroups", sep=",", file=fout)

    infolder = check_folder_path(infolder)
    for summary in glob.glob(infolder+'**/Statistics_Overall.csv', recursive=True):

        path_elem = summary.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]

        fin = open(summary, 'r')
        lines = fin.readlines()
        fin.close()

        summary_keywords = ["Number of genes", "Number of genes in orthogroups", "Number of unassigned genes",
                            "Number of orthogroups"]
        summary_values = []

        for su in summary_keywords:
            for line in lines:
                line_elem = line.strip().split('\t')
                if line_elem[0] == su:
                    summary_values.append(int(line_elem[-1]))

        summary_write = [str(s) for s in summary_values]
        print(genus, alignment, ",".join(summary_write), sep=",", file=fout)

    fout.close()


def generate_detailed_iadhore_summary(infolder, outfolder, gtf_folder,  genus_index=-3, alignment_index=-2):
    """
    Generate detailed i-ADHoRe summary
    :param gtf_folder: GTF folder
    :param infolder: Folder containing i-ADHoRe result
    :param outfolder: Output folder
    :param genus_index: Genus folder index from result
    :param alignment_index: Alignment folder index from result
    :return:
    """

    infolder = check_folder_path(infolder)
    for mult in glob.glob(infolder + '**/multiplicons.txt', recursive=True):

        path_elem = mult.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]
        subfolder = '/'.join([genus, alignment])

        outfolder = check_folder_path(outfolder, True)
        current_outfolder = outfolder + subfolder
        current_outfolder = check_folder_path(current_outfolder, True)
        outfile = current_outfolder + "iadhore_summary.csv"

        iadhore_result = os.path.dirname(mult)
        iadhore_dict = pi.iadhore_result_folder_to_dict(iadhore_result)
        iadhore_dict_with_location = cm.add_location_to_iadhore_synteny(iadhore_dict, gtf_folder)

        for key, value in iadhore_dict_with_location.items():
            for ke, val in value["segments"]:
                pass
