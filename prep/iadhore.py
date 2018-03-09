"""
i-ADHoRe Processing
"""

import glob
import os
import json
import yaml
import pickle

from misc.string import check_folder_path, is_int
from prep.gtf import gtf_to_dict_folder


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

            # get set of genes
            list_genes = set()
            fin = open(p, 'r')
            for line in fin.readlines():
                list_genes.add(line.strip()[:-1])
            fin.close()

            # check the difference
            exception_genes = list_genes.difference(family_genes)

            if len(exception_genes) < len(list_genes):

                sub_outfolder = outfolder + os.path.split(os.path.dirname(p))[-1]
                sub_outfolder = check_folder_path(sub_outfolder, True)

                outfile = sub_outfolder + os.path.basename(p)

                fout = open(outfile, 'w')

                # check infile again
                fin = open(p, 'r')
                for line in fin.readlines():
                    if line.strip()[:-1] not in exception_genes:
                        fout.write(line)
                fin.close()

                fout.close()


def create_iadhore_config(iadhore_genes_list, iadhore_family_file, iadhore_parameter_file, iadhore_result_folder,
                          outfile):
    """
    Create i-ADHoRe configuration file
    :param iadhore_genes_list: i-ADHoRe genes list
    :param iadhore_family_file: i-ADHoRe family file
    :param iadhore_parameter_file: i-ADHoRe parameter file
    :param iadhore_result_folder: i-ADHoRe result folder
    :param outfile: i-ADHoRe configuration file
    :return:
    """

    chromosome_dict = dict()
    for g in glob.glob(iadhore_genes_list + '**/*'):
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


def iadhore_result_file_to_dict(infile):
    """
    Convert i-ADHoRe result file (multiplicons.txt, segments.txt, list_elements.txt to dict)
    :param infile: i-ADHoRe result file
    :return: Python dictionary
    """

    results = dict()

    fin = open(infile, 'r')
    lines = fin.readlines()
    fin.close()

    attributes_keys = lines[0].strip().split("\t")
    for i in range(1, len(lines)):

        line = lines[i].strip()
        line_elem = line.split("\t")

        key = is_int(line_elem[0])

        attributes = dict()
        for j in range(1, len(line_elem)):
            attributes[attributes_keys[j]] = is_int(line_elem[j])

        results[key] = attributes

    return results


def merge_iadhore_file_dict(parent, child, parent_key_in_child, child_key_in_parent):
    """
    Merge 2 i-ADHoRe dictiory from i-ADHoRe output file
    :param parent: Parent dictionary
    :param child: Child dictionary
    :param parent_key_in_child: Parent key in child dictionary
    :param child_key_in_parent: Child key in parent dictionary
    :return: Merged dictionary
    """
    for key, value in parent.items():
        parent[key][child_key_in_parent] = dict()
        for k, v in child.items():
            if key == v[parent_key_in_child]:
                parent[key][child_key_in_parent][k] = v

    # delete the reference key
    for key in parent.keys():
        value = parent[key][child_key_in_parent]
        for k, v in value.items():
            del v[parent_key_in_child]

    return parent


def iadhore_result_folder_to_dict(infolder):
    """
    i-ADHoRe output folder to Python dictionary
    :param infolder: i-ADHoRe output folder
    :return: Python dictionary
    """

    infolder = check_folder_path(infolder)

    multiplicons_file = infolder + 'multiplicons.txt'
    segments_file = infolder + 'segments.txt'
    elements_file = infolder + 'list_elements.txt'

    # get dict
    multiplicons_dict = iadhore_result_file_to_dict(multiplicons_file)
    segments_dict = iadhore_result_file_to_dict(segments_file)
    elements_dict = iadhore_result_file_to_dict(elements_file)

    # merge segments dict with elements dict
    segments_dict = merge_iadhore_file_dict(segments_dict, elements_dict, "segment", "elements")

    # merge multiplicon dict with segment dict
    multiplicons_dict = merge_iadhore_file_dict(multiplicons_dict, segments_dict, "multiplicon", "segments")

    return multiplicons_dict


def get_complete_synteny_dict(input_dictionary):
    """
    Filter i-ADHoRe output folder dictionary by selecting only those which have orthologous in all species
    :param input_dictionary: i-ADHoRe output folder dictionary
    :return:
    """

    genomes = set()
    for value in input_dictionary.values():
        for v in value["segments"].values():
            genomes.add(v["genome"])
    number_of_genomes = len(genomes)

    complete_synteny_dict = dict()
    for key, value in input_dictionary.items():

        if len(value["segments"]) != number_of_genomes:
            continue

        genomes_in_segment = set()
        for v in value["segments"].values():
            genomes_in_segment.add(v["genome"])
        number_of_genomes_in_segment = len(genomes_in_segment)
        if number_of_genomes_in_segment != number_of_genomes:
            continue

        # get position intersection
        position_intersect = set()
        for v in value["segments"].values():

            position = set()
            for v0 in v["elements"].values():
                position.add(v0["position"])

            if not position_intersect:
                position_intersect = position
                continue

            position_intersect &= position

        new_segments = dict()
        for k, v in value["segments"].items():
            new_elements = dict()
            for k0, v0 in v["elements"].items():
                if v0["position"] in position_intersect:
                    new_elements[k0] = v0
            new_segments[k] = v
            new_segments[k]["elements"] = new_elements

        complete_synteny_dict[key] = value
        complete_synteny_dict[key]["segments"] = new_segments

    return complete_synteny_dict


def iadhore_result_to_serial(infolder, outfile, ftype="json", complete=False):
    """
    Convert i-ADHoRe output to data serialization file
    :param complete: Complete synteny, i.e., which has orthologous group in all species
    :param ftype: Type of the output file, e.g., JSON and YAML
    :param infolder: i-ADHoRe output folder
    :param outfile: Data serialization file
    :return:
    """

    serial_dict = iadhore_result_folder_to_dict(infolder)

    if complete:
        serial_dict = get_complete_synteny_dict(serial_dict)

    fout = open(outfile, 'w')

    if ftype=="json":
        json.dump(serial_dict, fout, indent=5)
    elif ftype=="yaml":
        yaml.dump(serial_dict, fout)
    elif ftype=="pickle":
        fout.close()
        fout = open(outfile, 'wb')
        pickle.dump(serial_dict, fout)

    fout.close()





