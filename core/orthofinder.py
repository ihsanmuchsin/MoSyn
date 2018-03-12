"""
OrthoFinder function
"""

import subprocess
import glob
import os

import prep.fasta as pf

from misc.string import check_folder_path


def run_orthofinder(infolder, outfolder, python2env, orthofinder_options=None):
    """
    Run OrthoFinder
    :param infolder: Input folder
    :param outfolder: Output folder
    :param python2env: Python 2 environment
    :param orthofinder_options: Other options for OrthoFinder
    :return:
    """
    outfolder = check_folder_path(outfolder, True)

    working_directory = outfolder + 'working_directory/'
    working_directory = check_folder_path(working_directory, True)

    cleaned_fasta = working_directory + 'cleaned_fasta/'
    cleaned_fasta = check_folder_path(cleaned_fasta, True)

    pf.clean_header_fasta_folder(infolder, cleaned_fasta)

    orthofinder_script = 'orthofinder.bash'
    fout = open(orthofinder_script, 'w')

    print('source', 'activate', python2env, file=fout)

    list_command = ['orthofinder', '-f', cleaned_fasta]
    if orthofinder_options:
        list_command += orthofinder_options.split(' ')
    print(' '.join(list_command), file=fout)

    fout.close()

    subprocess.call(['bash', orthofinder_script])
    subprocess.call(['rm', orthofinder_script])

    for p in glob.glob(cleaned_fasta + '*'):
        if os.path.isdir(p):

            p = check_folder_path(p)

            for r in glob.glob(p + '*'):
                subprocess.call(['mv', r, outfolder])

            subprocess.call(['rm', '-r', p])
