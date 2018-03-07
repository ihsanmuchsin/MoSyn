"""
The core of MoSyn pipeline
"""


from argparse import ArgumentParser
from sys import argv, exit

import core.orthofinder as co
import core.storm as cs
import core.iadhore as ci

from misc.string import check_folder_path


def orthofinder(args):

    infolder = check_folder_path(args.input)
    outfolder = check_folder_path(args.output)

    co.run_orthofinder(infolder, outfolder, python2env=args.python2, orthofinder_options=args.opt)


def storm(args):
    infolder = check_folder_path(args.input)
    pwm_folder = check_folder_path(args.pwm)
    outfolder = check_folder_path(args.output)

    storm_options = args.opt
    calculate_base_comp = args.bcomp

    cs.run_storm(infolder, pwm_folder, outfolder, storm_options, calculate_base_comp)


def iadhore(args):

    orthogroups_file = args.og
    gtf_folder = args.gtf
    id_conversion_folder = args.conv
    iadhore_parameter_file = args.param
    outfolder = args.output
    protein_column = args.pcol
    gene_column = args.gcol
    column_sep = args.csep

    ci.run_iadhore(orthogroups_file, gtf_folder, id_conversion_folder, iadhore_parameter_file, outfolder,
                protein_column, gene_column, column_sep)

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

p_st = subp.add_parser("storm", help="Run CREAD STORM")
p_st.add_argument("--input", metavar="<String>", help="Input folder", required=True)
p_st.add_argument("--pwm", metavar="<String>", help="PWM folder", required=True)
p_st.add_argument("--output", default="./storm_result/", metavar="<String>",
                  help="Output folder")
p_st.add_argument("--opt", default=None, metavar="<String>",
                  help="Options for CREAD STORM. "
                       "Must be put inside a single quote, for example '-t 17'. "
                       "Please refer to the CREAD STORM Manual.")
p_st.add_argument("--bcomp", default=False, action='store_true', help="Calculate base composition")
p_st.set_defaults(func=storm)

p_ia = subp.add_parser("i-adhore", help="Run i-ADHoRe")
p_ia.add_argument("--og", metavar="<String>", help="Orthogroups file", required=True)
p_ia.add_argument("--gtf", metavar="<String>", help="GTF folder", required=True)
p_ia.add_argument("--conv", metavar="<String>", help="ID conversion folder", required=True)
p_ia.add_argument("--param", metavar="<String>", help="i-ADHoRe parameter file", required=True)
p_ia.add_argument("--output", default="./iadhore_result/", metavar="<String>", help="Output folder")
p_ia.add_argument("--pcol", default=0, metavar="<Integer>", help="Protein column")
p_ia.add_argument("--gcol", default=1, metavar="<Integer>", help="Gene column")
p_ia.add_argument("--csep", default="\t", metavar="<String>", help="Column separator")
p_ia.set_defaults(func=iadhore)


if len(argv) == 1:
    p.print_help()
    exit(0)

if len(argv) == 2:
    if argv[1] == 'orthofinder':
        p_of.print_help()
    elif argv[1] == 'storm':
        p_st.print_help()
    elif argv[1] == 'i-adhore':
        p_ia.print_help()
    exit(0)

args = p.parse_args()

args.func(args)
