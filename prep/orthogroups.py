"""
OrthoFinder/OrthoMCL Orthogroups File Processing
"""

from collections import Counter

import prep.idconversion as pcv


def orthogroups_file_to_dict(infile):
    """
    Get orthogroups dictionary from orthogroups file
    :param infile: Orthogroups file
    :return: Python dictionary of the orthogroups
    """

    orthogroups_dict = dict()

    fin = open(infile, 'r')

    for line in fin.readlines():
        key, value = line.strip().split(":")
        value = value.strip().split(" ")
        orthogroups_dict[key] = value

    fin.close()

    return orthogroups_dict


def orthogroups_protein_to_gene(orthogroups_file, id_conversion_folder, outfile, protein_column=0, gene_column=1, column_sep="\t", non_redundant=True):
    """
    Convert protein ID orthogroups to gene ID orthogroups
    :param non_redundant: Remove redundancies
    :param orthogroups_file: protein ID orthogroups file
    :param id_conversion_folder: ID conversion file
    :param outfile: gene ID orthogroups file
    :param protein_column: Column containing protein ID
    :param gene_column: Column containing gene ID
    :param column_sep: Separator
    :return:
    """

    orthogroups_dict = orthogroups_file_to_dict(orthogroups_file)
    id_conversion_dict = pcv.id_conversion_folder_to_dict(id_conversion_folder, protein_column, gene_column, column_sep)

    converted_orthogroups_dict = dict()
    for key, values in orthogroups_dict.items():
        converted_values = set()
        for value in values:
            converted_values.add(id_conversion_dict[value])
        converted_values = list(converted_values)
        converted_orthogroups_dict[key] = converted_values

    # Non redundant removes value redundancies
    if non_redundant:
        # Get all value
        all_values = []
        for values in converted_orthogroups_dict.values():
            all_values.extend(values)

        # Count each value occurence and get the list of values which occur more than 1
        value_counter = Counter(all_values)
        duplicate_values = [k for k, v in value_counter.items() if v > 1]

        # Get duplicate value to keys
        duplicate_dict = dict()
        for key, values in converted_orthogroups_dict.items():
            for value in values:
                if value in duplicate_values:
                    if value not in duplicate_dict.keys():
                        duplicate_dict[value] = []
                    duplicate_dict[value].append(key)

        # remove duplicate value, only save one
        for key, values in duplicate_dict.items():
            for i in range(1, len(values)):
                converted_orthogroups_dict[values[i]].remove(key)

        # remove empty value from final dictionary
        empty_values = []
        for key, values in converted_orthogroups_dict.items():
            if not values:
                empty_values.append(key)

        for key in empty_values:
            del converted_orthogroups_dict[key]

    # Now save the dictionary to file
    fout = open(outfile, 'w')
    for key, values in converted_orthogroups_dict.items():
        print (key+":", " ".join(values), file=fout)
    fout.close()


def orthogroups_to_iadhore_family_dict(infile):
    """
    Convert Orthogroups to i-ADHoRe family dictionary
    :param infile: Orthogroups file
    :return: i-ADHoRe family dictionary
    """

    iadhore_family_dict = dict()

    orthogroups_dict = orthogroups_file_to_dict(infile)
    for key, values in orthogroups_dict.items():
        for value in values:
            iadhore_family_dict[value] = key

    return iadhore_family_dict


def orthogroups_to_iadhore_family_file(infile, outfile):
    """
    Convert Orthogroups to i-ADHoRe family file
    :param infile: Orthogroups file
    :param outfile: i-ADHoRe family file
    :return:
    """

    fout = open(outfile, 'w')

    iadhore_family_dict = orthogroups_to_iadhore_family_dict(infile)
    for key, value in iadhore_family_dict.items():
        print(key, value, sep="\t", file=fout)

    fout.close()