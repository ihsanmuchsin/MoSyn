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


def iadhore_list_family_filtering(iadhore_list, iadhore_family_file, outfolder):
    """
    Select only genes that are in the family file
    :param outfolder: Output folder
    :param iadhore_list: i-ADHoRe list
    :param iadhore_family_file: i-ADHoRe family file
    :return:
    """

    iadhore_list = check_folder_path(iadhore_list)
    outfolder = check_folder_path(outfolder, True)

    family_genes = set(iadhore_family_to_dict(iadhore_family_file).keys())

    for p in glob.glob(iadhore_list + '**/*'):
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