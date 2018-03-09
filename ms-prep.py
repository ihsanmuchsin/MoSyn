"""
Data pre processing for MoSyn pipeline
"""

from argparse import ArgumentParser
from sys import argv, exit

import prep.fasta as pf
import prep.gtf as pg
import prep.orthofinder as pof
import prep.iadhore as pi
import prep.storm as ps


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

    pof.orthogroups_protein_to_gene(orthogroup_file, id_conversion_folder, outfile, protein_column, gene_column, column_sep, non_redundant)


def og2fam(args):
    infile = args.input
    outfile = args.output
    pof.orthogroups_to_iadhore_family_file(infile, outfile)


def gtf2list(args):
    infolder = args.input
    outfolder = args.output
    pg.gtf_to_iadhore_list_folder(infolder, outfolder)


def listfilter(args):
    infolder = args.input
    family_file = args.fam
    outfolder = args.output

    pi.iadhore_list_family_filtering(infolder, family_file, outfolder)


def createiaconfig(args):
    iadhore_genes_list = args.glist
    iadhore_family_file = args.fam
    iadhore_parameter_file = args.param
    iadhore_result_folder = args.res
    outfile = args.output

    pi.create_iadhore_config(iadhore_genes_list, iadhore_family_file, iadhore_parameter_file, iadhore_result_folder, outfile)


def iadhore2serial(args):

    infolder = args.input
    outfile = args.output

    pi.iadhore_result_to_serial(infolder, outfile, args.ftype, args.complete)

def transfac2gtf(args):

    infolder = args.input
    outfolder = args.output
    ps.transfac_to_gtf_folder(infolder, outfolder, args.mid, args.idx, args.name)


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
p_optg.add_argument('--csep', metavar='<String>', help='Column separator', default="\t")
p_optg.add_argument('--nr', help='Non redundant values', default=False, action='store_true')
p_optg.set_defaults(func=ogprot2gene)

p_otf = subp.add_parser('og2fam', help='Convert Orthogroups file to i-ADHoRe family')
p_otf.add_argument('--input', metavar='<String>', help='Path to the input file', required=True)
p_otf.add_argument('--output', metavar='<String>', help='Path to the output file', required=True)
p_otf.set_defaults(func=og2fam)

p_gtl = subp.add_parser('gtf2list', help='Convert GTF files to i-ADHoRe genes list')
p_gtl.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_gtl.add_argument('--output', metavar='<String>', help='Path to the output folder', required=True)
p_gtl.set_defaults(func=gtf2list)

p_lf = subp.add_parser('listfilter', help='Filter i-ADHoRe genes list by choosing only genes in orthogroups')
p_lf.add_argument('--input', metavar='<String>', help='Path to the i-ADHoRe genes list', required=True)
p_lf.add_argument('--fam', metavar='<String>', help='Path to the i-ADHoRe family file', required=True)
p_lf.add_argument('--output', metavar='<String>', help='Path to the output folder', required=True)
p_lf.set_defaults(func=listfilter)

p_cc = subp.add_parser('createiaconfig', help='Create i-ADHoRe configuration file')
p_cc.add_argument('--glist', metavar='<String>', help='Path to the i-ADHoRe genes list', required=True)
p_cc.add_argument('--fam', metavar='<String>', help='Path to the i-ADHoRe family file', required=True)
p_cc.add_argument('--param', metavar='<String>', help='Path to the i-ADHoRe parameter file', required=True)
p_cc.add_argument('--res', metavar='<String>', help='Path to the i-ADHoRe output folder', required=True)
p_cc.add_argument('--output', metavar='<String>', help='Path to the output file', required=True)
p_cc.set_defaults(func=createiaconfig)

p_iti = subp.add_parser('ia2serial', help='Convert i-ADHoRe output to JSON/YAML file')
p_iti.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_iti.add_argument('--output', metavar='<String>', help='Path to the output file', required=True)
p_iti.add_argument('--ftype', metavar='<String>', help='Type of the output file', default="json")
p_iti.add_argument('--complete', help='Complete synteny', default=False, action='store_true')
p_iti.set_defaults(func=iadhore2serial)

p_stg = subp.add_parser('transfac2gtf', help='Convert STORM output / TRANSFAC file format to GTF')
p_stg.add_argument('--input', metavar='<String>', help='Path to the input folder', required=True)
p_stg.add_argument('--output', metavar='<String>', help='Path to the output folder', default="./STORM_to_GTF")
p_stg.add_argument('--mid', metavar='<String>', help='Binding sites ID', default="BS")
p_stg.add_argument('--idx', metavar='<String>', help='The starting number of the ID index', default=0, type=int)
p_stg.add_argument('--name', metavar='<String>', help='The Motif name', default="MOTIF")
p_stg.set_defaults(func=transfac2gtf)


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
    elif argv[1] == 'og2fam':
        p_otf.print_help()
    elif argv[1] == 'gtf2list':
        p_gtl.print_help()
    elif argv[1] == 'listfilter':
        p_lf.print_help()
    elif argv[1] == 'createiaconfig':
        p_cc.print_help()
    elif argv[1] == 'ia2serial':
        p_iti.print_help()
    elif argv[1] == 'transfac2gtf':
        p_stg.print_help()
    exit(0)

args = p.parse_args()
args.func(args)
