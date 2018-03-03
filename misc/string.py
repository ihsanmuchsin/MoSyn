"""
Miscellaneous methods for MoSyn pipeline
"""

import os


def check_folder_path(folder, create_dir=False):
    """

    :param folder: Check if the last character of the path is /. This also create a new folder if it does not exist
    :return: Corrected path of the folder
    """

    if create_dir:
        if not os.path.exists(folder):
            os.makedirs(folder)

    if folder[-1] != "/":
        folder += "/"

    return folder