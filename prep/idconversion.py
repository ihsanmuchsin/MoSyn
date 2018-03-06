"""
ID Conversion File Processing
"""

import glob

from misc.string import check_folder_path


def id_conversion_file_to_dict(infile, protein_column=0, gene_column=1, column_sep="\t"):
    """
    Convert Protein ID -> Gene ID conversion file to dictionary
    :param infile: ID conversion file
    :param protein_column: Column containing protein ID
    :param gene_column: Column containing gene ID
    :param column_sep: separator
    :return:
    """

    id_conversion_dict = dict()

    fin = open(infile, 'r')

    for line in fin.readlines():
        if not line.startswith("#"):

            line_elem = line.strip().split(column_sep)
            protein_id = line_elem[protein_column]
            gene_id = line_elem[gene_column]
            id_conversion_dict[protein_id] = gene_id

    fin.close()

    return id_conversion_dict


def id_conversion_folder_to_dict(infolder, protein_column=0, gene_column=1, column_sep="\t"):

    super_id_conversion_dict = dict()

    infolder = check_folder_path(infolder)

    for file in glob.glob(infolder + '*'):
        id_conversion_dict = id_conversion_file_to_dict(file, protein_column, gene_column, column_sep)
        for key, value in id_conversion_dict.items():
            super_id_conversion_dict[key] = value

    return super_id_conversion_dict