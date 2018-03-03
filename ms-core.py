"""
The core of MoSyn pipeline
"""


from argparse import ArgumentParser
from sys import argv, exit

from core.orthofinder import *
from misc.string import *


def orthofinder(args):

    infolder = check_folder_path(args.input)
    outfolder = check_folder_path(args.output)

    run_orthofinder(infolder, outfolder, python2env=args.python2, orthofinder_options=args.opt)


# create parser
p = ArgumentParser(prog='ms-core', description='MoSyn Pipeline')

subp = p.add_subparsers()

p_of = subp.add_parser("orthofinder", help="Run OrthoFinder")
p_of.add_argument("--input", metavar="<String>", help="Input folder", required=True)
p_of.add_argument("--output", default="./orthofinder_result/", metavar="<String>",
                  help="Output folder")
p_of.add_argument("--python2", metavar='<String>', help="Python 2 environment", required=True)
p_of.add_argument("--opt", default=None, metavar="<String>",
                  help="Options for OrthoFinder. "
                       "Must be put inside a single quote, for example '-S diamond'. "
                       "Please refer to the OrthoFinder Manual.")
p_of.set_defaults(func=orthofinder)


if len(argv) == 1:
    p.print_help()
    exit(0)

if len(argv) == 2:
    if argv[1] == 'orthofinder':
        p_of.print_help()
    exit(0)

args = p.parse_args()

args.func(args)