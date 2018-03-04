"""
Data pre processing for MoSyn pipeline
"""

from argparse import ArgumentParser
from sys import argv, exit

import preprocess.fasta as pf
import preprocess.gtf as pg


def fasta2json(args):
    infolder = args.input
    outfolder = args.output
    pf.fasta_to_json_folder(infolder, outfolder, args.type)


def cleanheader(args):
    infolder = args.input
    outfolder = args.output
    pf.clean_header_fasta_folder(infolder, outfolder)


def gtf2json(args):
    infolder = args.input
    outfolder = args.output
    pg.gtf_to_json_folder(infolder, outfolder)


p = ArgumentParser(prog='ms-prep', description='Data preprocessing for MoSyn pipeline')

subp = p.add_subparsers()

p_ftj = subp.add_parser('fasta2json', help='Convert FASTA to JSON')
p_ftj.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_ftj.add_argument('--output', metavar='<String>', help='Path to the output folder', required=True)
p_ftj.add_argument('--type', metavar='<String>', help='Type of the input', required=True)
p_ftj.set_defaults(func=fasta2json)

p_cfh = subp.add_parser('cleanheader', help='Clean FASTA header')
p_cfh.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_cfh.add_argument('--output', metavar='<String>', help='Path to the output folder', required=True)
p_cfh.set_defaults(func=cleanheader)

p_gtj = subp.add_parser('gtf2json', help='Convert GTF to JSON')
p_gtj.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_gtj.add_argument('--output', metavar='<String>', help='Path to the output folder', required=True)
p_gtj.set_defaults(func=gtf2json)

if len(argv) == 1:
    p.print_help()
    exit(0)

if len(argv) == 2:
    if argv[1] == 'fasta2json':
        p_ftj.print_help()
    elif argv[1] == 'cleanheader':
        p_cfh.print_help()
    exit(0)

args = p.parse_args()
args.func(args)
