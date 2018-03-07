"""
i-ADHoRe function
"""

import subprocess

from prep.orthofinder import orthogroups_protein_to_gene, orthogroups_to_iadhore_family_file
from prep.gtf import gtf_to_iadhore_list_folder
from prep.iadhore import iadhore_list_family_filtering, create_iadhore_config
from misc.string import check_folder_path


def run_iadhore(orthogroups_file, gtf_folder, id_conversion_folder, iadhore_parameter_file, outfolder,
                protein_column=0, gene_column=1, column_sep="\t"):
    """
    Run i-ADHoRe
    :param orthogroups_file: Orthogroups file
    :param gtf_folder: GTF folder
    :param id_conversion_folder: ID conversion folder
    :param iadhore_parameter_file: i-ADHoRe parameter file
    :param outfolder: i-ADHoRe output folder
    :param protein_column: Protein column
    :param gene_column: Gene column
    :param column_sep: Column separator
    :return:
    """

    outfolder = check_folder_path(outfolder, True)

    working_directory = outfolder + 'working_directory/'
    working_directory = check_folder_path(working_directory, True)

    orthogroups_genes = working_directory + 'Orthogroups.txt'
    orthogroups_protein_to_gene(orthogroups_file, id_conversion_folder, orthogroups_genes,
                                protein_column, gene_column, column_sep)

    iadhore_family_file = working_directory + 'iadhore_family.tsv'
    orthogroups_to_iadhore_family_file(orthogroups_genes, iadhore_family_file)

    all_genes_list = working_directory + 'temporary_genes_list/'
    all_genes_list = check_folder_path(all_genes_list, True)
    gtf_to_iadhore_list_folder(gtf_folder, all_genes_list)

    filtered_genes_list = working_directory + 'genes_list'
    filtered_genes_list = check_folder_path(filtered_genes_list, True)
    iadhore_list_family_filtering(all_genes_list, iadhore_family_file, filtered_genes_list)

    # remove temporary genes list
    subprocess.call(['rm', '-r', all_genes_list])

    iadhore_config = working_directory + 'iadhore_config.ini'
    create_iadhore_config(filtered_genes_list, iadhore_family_file, iadhore_parameter_file, outfolder, iadhore_config)

    # run i-ADHoRe
    subprocess.call(['i-adhore', iadhore_config])
