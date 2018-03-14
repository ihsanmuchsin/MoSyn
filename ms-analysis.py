"""
The analysis of MoSyn output
"""


from argparse import ArgumentParser
from sys import argv, exit

import yaml

import analysis.loop as alo
import analysis.summary as asm

from misc.string import check_folder_path


def idloop(args):

    outdir = args.output
    outdir = check_folder_path(outdir, True)
    outfile = outdir + "loops.csv"
    outfile0 = outdir + "loops.yaml"

    infile = args.input
    mosyn_loops = alo.identify_loops_in_synteny(infile, args.min, args.max)
    alo.print_loops_to_csv(mosyn_loops, outfile)
    with open(outfile0, 'w') as stream:
        yaml.dump(mosyn_loops, stream)

def summarize_orthofinder(args):

    storm_result = args.input
    outfile = args.output

    asm.generate_orthofinder_summary(storm_result, outfile, args.gidx, args.aidx)


p = ArgumentParser(prog='ms-analysis', description='Data analysis for MoSyn pipeline')

subp = p.add_subparsers()

p_ilo = subp.add_parser('idloop', help='Identify loops from synteny.yaml file')
p_ilo.add_argument('--input', metavar='<String>', help='Path to the input file', required=True)
p_ilo.add_argument('--output', metavar='<String>', help='Path to the output folder', default="./loops/")
p_ilo.add_argument('--min', metavar='<Integer>', help='The minimum length of the loops', default=80000, type=int)
p_ilo.add_argument('--max', metavar='<Integer>', help='The maximum length of the loops', default=800000, type=int)
p_ilo.set_defaults(func=idloop)

p_suo = subp.add_parser('summarize-ortho', help='Summarize OrthoFinder result')
p_suo.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_suo.add_argument('--output', metavar='<String>', help='Path to the output file', default="orthofinder_summary.csv")
p_suo.add_argument('--gidx', metavar='<Integer>', help='The minimum length of the loops', default=-3, type=int)
p_suo.add_argument('--aidx', metavar='<Integer>', help='The maximum length of the loops', default=-2, type=int)
p_suo.set_defaults(func=summarize_orthofinder)


if len(argv) == 1:
    p.print_help()
    exit(0)

if len(argv) == 2:
    if argv[1] == 'idloop':
        p_ilo.print_help()
    elif argv[1] == 'summarize-ortho':
        p_suo.print_help()
    exit(0)

args = p.parse_args()
args.func(args)
