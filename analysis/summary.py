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

    print("Genus", "Alignment", "Number_of_Genes", "Number_of_Genes_in_Orthogroups", "Number_of_Unassigned_Genes",
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


def generate_detailed_iadhore_summary(infolder, outfile, material_folder, genus_index=-3, alignment_index=-2):
    """
    Generate detailed i-ADHoRe summary.
    Result structure Genus -> Alignment -> Result
    GTF structure Material -> Genus -> GTF
    :param material_folder: Folder containing GTF folder
    :param infolder: Folder containing i-ADHoRe result
    :param outfile: Output folder
    :param genus_index: Genus folder index from result
    :param alignment_index: Alignment folder index from result
    :return:
    """

    material_folder = check_folder_path(material_folder)

    fout = open(outfile, 'w')

    print("Genus", "Alignment", "Synteny_ID", "Number_of_Segments", "Number_of_Genes", "Total_Length",
          sep=",", file=fout)

    infolder = check_folder_path(infolder)
    for mult in glob.glob(infolder + '**/multiplicons.txt', recursive=True):

        path_elem = mult.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]

        gtf_folder = material_folder + "/".join([genus, "GTF"])
        gtf_folder = check_folder_path(gtf_folder)

        iadhore_result = os.path.dirname(mult)
        iadhore_dict = pi.iadhore_result_folder_to_dict(iadhore_result)
        iadhore_dict_with_location = cm.add_location_to_iadhore_synteny(iadhore_dict, gtf_folder)

        for key, value in iadhore_dict_with_location.items():

            number_of_segments = 0
            number_of_genes = 0
            total_length = 0

            for ke, val in value["segments"].items():
                number_of_segments += 1
                length = abs(val["start"] - val["end"])
                total_length += length
                number_of_genes += len(val["elements"])

            print(genus, alignment, key, number_of_segments, number_of_genes, total_length, sep=",", file=fout)

    fout.close()


def generate_short_iadhore_summary(infolder, outfile, material_folder, genus_index=-3, alignment_index=-2):
    """
    Generate detailed i-ADHoRe summary.
    Result structure Genus -> Alignment -> Result
    GTF structure Material -> Genus -> GTF
    :param material_folder: Folder containing GTF folder
    :param infolder: Folder containing i-ADHoRe result
    :param outfile: Output folder
    :param genus_index: Genus folder index from result
    :param alignment_index: Alignment folder index from result
    :return:
    """

    material_folder = check_folder_path(material_folder)

    fout = open(outfile, 'w')

    print("Genus", "Alignment", "Number_of_Synteny", "Number_of_Genes_in_Synteny",
          "Average_Number_of_Genes_per_Synteny", "Average_Length_per_Synteny",
          sep=",", file=fout)

    infolder = check_folder_path(infolder)
    for mult in glob.glob(infolder + '**/multiplicons.txt', recursive=True):

        path_elem = mult.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]

        gtf_folder = material_folder + "/".join([genus, "GTF"])
        gtf_folder = check_folder_path(gtf_folder)

        iadhore_result = os.path.dirname(mult)
        iadhore_dict = pi.iadhore_result_folder_to_dict(iadhore_result)
        iadhore_dict_with_location = cm.add_location_to_iadhore_synteny(iadhore_dict, gtf_folder)

        synteny_length = 0
        synteny_genes = 0
        num_of_mult = 0

        set_of_genes = set()

        for key, value in iadhore_dict_with_location.items():

            number_of_segments = 0
            mult_genes = 0
            mult_length = 0

            for ke, val in value["segments"].items():
                number_of_segments += 1
                length = abs(val["start"] - val["end"])
                mult_length += length
                mult_genes += len(val["elements"])

                for v in val["elements"].values():
                    set_of_genes.add(v["gene"])

            mult_avg_len = mult_length / number_of_segments
            mult_avg_genes = mult_genes / number_of_segments

            synteny_length += mult_avg_len
            synteny_genes += mult_avg_genes
            num_of_mult += 1

        nr_genes = len(set_of_genes)
        synteny_avg_length = synteny_length / num_of_mult
        synteny_avg_genes = synteny_genes / num_of_mult

        print(genus, alignment, num_of_mult, nr_genes, synteny_avg_genes, synteny_avg_length, sep=",", file=fout)

    fout.close()

