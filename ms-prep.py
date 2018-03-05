"""
Data pre processing for MoSyn pipeline
"""

from argparse import ArgumentParser
from sys import argv, exit

import prep.fasta as pf
import prep.gtf as pg
import prep.orthogroups as pog


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


def fasta2idconv(args):
    infolder = args.input
    outfolder = args.output
    pf.fasta_to_conversion_folder(infolder, outfolder, args.type)


def ogprot2gene(args):
    orthogroup_file = args.og
    id_conversion_folder = args.conv
    outfile = args.output

    protein_column = args.pcol
    gene_column = args.gcol
    column_sep = args.csep

    non_redundant = args.nr

    pog.orthogroups_protein_to_gene(orthogroup_file, id_conversion_folder, outfile, protein_column, gene_column, column_sep, non_redundant)


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

p_fti = subp.add_parser('fasta2idconv', help='Convert FASTA to ID Conversion file')
p_fti.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_fti.add_argument('--output', metavar='<String>', help='Path to the output folder', required=True)
p_fti.add_argument('--type', metavar='<String>', help='Type of the input', required=True)
p_fti.set_defaults(func=fasta2idconv)

p_optg = subp.add_parser('ogprot2gene', help='Convert Orthogroups Protein to Gene')
p_optg.add_argument('--og', metavar='<String>', help='Orthogroups protein file', required=True)
p_optg.add_argument('--conv', metavar='<String>', help='Protein ID to Gene ID conversion folder', required=True)
p_optg.add_argument('--output', metavar='<String>', help='Path to the output file', required=True)
p_optg.add_argument('--pcol', metavar='<Integer>', help='Protein column', type=int, default=0)
p_optg.add_argument('--gcol', metavar='<Integer>', help='Gene column', type=int, default=1)
p_optg.add_argument('--csep', metavar='<Integer>', help='Column separator', default="\t")
p_optg.add_argument('--nr', help='Non redundant values', default=False, action='store_true')
p_optg.set_defaults(func=ogprot2gene)

if len(argv) == 1:
    p.print_help()
    exit(0)

if len(argv) == 2:
    if argv[1] == 'fasta2json':
        p_ftj.print_help()
    elif argv[1] == 'cleanheader':
        p_cfh.print_help()
    elif argv[1] == 'gtf2json':
        p_gtj.print_help()
    elif argv[1] == 'fasta2idconv':
        p_fti.print_help()
    elif argv[1] == 'ogprot2gene':
        p_optg.print_help()
    exit(0)

args = p.parse_args()
args.func(args)
