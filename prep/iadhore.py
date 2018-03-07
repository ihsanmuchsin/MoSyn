"""
i-ADHoRe Processing
"""

import glob
import os

from misc.string import check_folder_path


def iadhore_family_to_dict(infile):
    """
    Convert i-ADHoRe family to Python dictionary
    :param infile: i-ADHoRe family file
    :return: Python dictionary Gene ID -> Family ID
    """

    family_dict = dict()

    fin = open(infile, 'r')
    for line in fin.readlines():
        key, value = line.strip().split("\t")
        family_dict[key] = value

    return family_dict


def iadhore_list_family_filtering(iadhore_genes_list, iadhore_family_file, outfolder):
    """
    Select only genes that are in the family file
    :param outfolder: Output folder
    :param iadhore_genes_list: i-ADHoRe list
    :param iadhore_family_file: i-ADHoRe family file
    :return:
    """

    iadhore_genes_list = check_folder_path(iadhore_genes_list)
    outfolder = check_folder_path(outfolder, True)

    family_genes = set(iadhore_family_to_dict(iadhore_family_file).keys())

    for p in glob.glob(iadhore_genes_list + '**/*'):
        if os.path.isfile(p):

            #get set of genes
            list_genes = set()
            fin = open(p, 'r')
            for line in fin.readlines():
                list_genes.add(line.strip()[:-1])
            fin.close()

            #check the difference
            exception_genes = list_genes.difference(family_genes)

            if len(exception_genes) < len(list_genes):

                sub_outfolder = outfolder + os.path.split(os.path.dirname(p))[-1]
                sub_outfolder = check_folder_path(sub_outfolder, True)

                outfile = sub_outfolder + os.path.basename(p)

                fout = open(outfile, 'w')

                #check infile again
                fin = open(p, 'r')
                for line in fin.readlines():
                    if line.strip()[:-1] not in exception_genes:
                        fout.write(line)
                fin.close()

                fout.close()


def create_iadhore_config(iadhore_genes_list, iadhore_family_file, iadhore_parameter_file, iadhore_result_folder, outfile):

    chromosome_dict = dict()
    for g in glob.glob(iadhore_genes_list+ '**/*'):
        if os.path.isfile(g):
            species = g.split('/')[-2]
            if species not in chromosome_dict.keys():
                chromosome_dict[species] = []
                chromosome_dict[species].append(g)
            else:
                chromosome_dict[species].append(g)

    fout = open(outfile, 'w')
    for key in chromosome_dict.keys():
        fout.write('genome=' + key + '\n')
        for val in chromosome_dict[key]:
            chromosome = os.path.basename(val).split('.')[0]
            fout.write(chromosome + ' ' + val + '\n')
        fout.write('\n')

    iadhore_result_folder = check_folder_path(iadhore_result_folder, True)

    fout.write('output_path=' + iadhore_result_folder + '\n')
    fout.write('blast_table=' + iadhore_family_file + '\n')
    fout.write('table_type=family\n')

    with open(iadhore_parameter_file, 'r') as fin:
        for line in fin.readlines():
            fout.write(line)

    fout.close()
