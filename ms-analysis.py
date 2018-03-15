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


def summarize(args):

    infolder = args.input

    if args.type == "orthofinder":
        asm.generate_orthofinder_summary(infolder, args.gidx, args.aidx)
    elif args.type == "mosyn":
        asm.generate_mosyn_summary(infolder, args.gidx, args.aidx, args.sidx, args.pidx)
    elif args.type == "loop":
        asm.generate_loop_summary(infolder, args.output, args.gidx, args.aidx, args.sidx, args.pidx)


p = ArgumentParser(prog='ms-analysis', description='Data analysis for MoSyn pipeline')

subp = p.add_subparsers()

p_ilo = subp.add_parser('idloop', help='Identify loops from synteny.yaml file')
p_ilo.add_argument('--input', metavar='<String>', help='Path to the input file', required=True)
p_ilo.add_argument('--output', metavar='<String>', help='Path to the output folder', default="./loops/")
p_ilo.add_argument('--min', metavar='<Integer>', help='The minimum length of the loops', default=80000, type=int)
p_ilo.add_argument('--max', metavar='<Integer>', help='The maximum length of the loops', default=800000, type=int)
p_ilo.set_defaults(func=idloop)

p_suo = subp.add_parser('summarize', help='Summarize result')
p_suo.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_suo.add_argument('--type', metavar='<String>', help='Type of the result to summarize, e.g., orthofinder',
                   required=True)
p_suo.add_argument('--output', metavar='<String>', default="./summary/",
                   help='Path to the output folder. Required for Loops summary')
p_suo.add_argument('--material', metavar='<String>', help='Path to the material folder. Required for i-ADHoRe summary')
p_suo.add_argument('--gidx', metavar='<Integer>', help='The genus folder relative index to result', type=int)
p_suo.add_argument('--aidx', metavar='<Integer>', help='The alignment folder relative index to result', type=int)
p_suo.add_argument('--sidx', metavar='<Integer>', help='The score folder relative index to result', type=int)
p_suo.add_argument('--pidx', metavar='<Integer>', help='The pwm folder relative index to result', type=int)
p_suo.set_defaults(func=summarize)


if len(argv) == 1:
    p.print_help()
    exit(0)

if len(argv) == 2:
    if argv[1] == 'idloop':
        p_ilo.print_help()
    elif argv[1] == 'summarize':
        p_suo.print_help()
    exit(0)

args = p.parse_args()
args.func(args)
