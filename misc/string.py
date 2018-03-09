"""
Miscellaneous methods for MoSyn pipeline
"""

import os


def check_folder_path(folder, create_dir=False):
    """

    :param create_dir: Create new directory
    :param folder: Check if the last character of the path is /. This also create a new folder if it does not exist
    :return: Corrected path of the folder
    """

    if create_dir:
        if not os.path.exists(folder):
            os.makedirs(folder)

    if folder[-1] != "/":
        folder += "/"

    return folder


def is_float(float_string):
    """
    Check if the string can be converted to float. Return float if TRUE
    :param float_string: String
    :return: Float or String
    """
    try:
        float_string = float(float_string)
    except ValueError:
        return float_string
    return float_string


def is_int(int_string):
    """
    Check if the string can be converted to int. Return int if TRUE
    :param int_string: String
    :return: Integer or String
    """
    try:
        int_string = int(int_string)
    except ValueError:
        return int_string
    return int_string
