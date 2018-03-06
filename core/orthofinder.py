"""
OrthoFinder function
"""

import subprocess
import glob
import os

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
    orthofinder_script = 'orthofinder.bash'
    fout = open(orthofinder_script, 'w')

    print('source', 'activate', python2env, file=fout)

    list_command = ['orthofinder', '-f', infolder]
    if orthofinder_options != None:
        list_command += orthofinder_options.split(' ')
    print(' '.join(list_command), file=fout)

    fout.close()

    subprocess.call(['bash', orthofinder_script])
    subprocess.call(['rm', orthofinder_script])

    outfolder = check_folder_path(outfolder, True)

    for p in glob.glob(infolder + '*'):
        if os.path.isdir(p):

            p = check_folder_path(p)

            for r in glob.glob(p + '*'):
                subprocess.call(['mv', r, outfolder])

            subprocess.call(['rm', '-r', p])
