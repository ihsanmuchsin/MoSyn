"""
The core of MoSyn pipeline
"""


from argparse import ArgumentParser
from sys import argv, exit

import glob
import os

import core.orthofinder as co
import core.storm as cs
import core.iadhore as ci
import core.mosyn as cm

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


def mosyn(args):
    iadhore_output_folder = args.iadhore
    storm_output_folder = args.storm
    gtf_folder = args.gtf
    outfolder = args.output
    window = args.window
    complete = args.complete
    binding_site_id = args.mid
    id_start_index = args.idx
    motif_name = args.name

    cm.run_mosyn(iadhore_output_folder, storm_output_folder, gtf_folder, outfolder,
                 window, complete, binding_site_id, id_start_index, motif_name)


def runall(args):

    outfolder = args.output
    outfolder = check_folder_path(outfolder, True)

    working_directory = outfolder + "working_directory/"
    working_directory = check_folder_path(working_directory, True)

    proteome = args.proteome
    orthofinder_output = working_directory + "orthofinder_output/"
    orthofinder_output = check_folder_path(orthofinder_output, True)

    co.run_orthofinder(proteome, orthofinder_output,
                       python2env=args.python2, orthofinder_options=args.ofopt)

    gtf_folder = args.gtf
    id_conversion_folder = args.idconv
    iadhore_parameter_file = args.param
    orthogroups_file = orthofinder_output + "Orthogroups.txt"
    iadhore_output = working_directory + "iadhore_output/"
    iadhore_output = check_folder_path(iadhore_output, True)

    ci.run_iadhore(orthogroups_file, gtf_folder, id_conversion_folder, iadhore_parameter_file, iadhore_output,
                   args.pcol, args.gcol, args.csep)

    genome_folder = args.genome
    pwm_folder = check_folder_path(args.pwm)

    storm_output = working_directory + "storm_output/"
    storm_output = check_folder_path(storm_output, True)

    cs.run_storm(genome_folder, pwm_folder, storm_output,
                 storm_options=args.stopt, calculate_base_comp=args.bcomp)

    for pwm in glob.glob(pwm_folder + '*'):
        sub_storm = storm_output + os.path.basename(pwm).split('.')[0]
        sub_storm = check_folder_path(sub_storm)

        sub_out = outfolder + os.path.basename(pwm).split('.')[0]
        sub_out = check_folder_path(sub_out, True)

        cm.run_mosyn(iadhore_output, sub_storm, gtf_folder, sub_out,
                     args.window, args.complete, args.mid, args.idx, args.name)


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

p_ms = subp.add_parser('mosyn', help='Run MoSyn')
p_ms.add_argument('--iadhore', metavar='<String>', help='Path to the iadhore result folder', required=True)
p_ms.add_argument('--storm', metavar='<String>', help='Path to the storm result folder', required=True)
p_ms.add_argument('--gtf', metavar='<String>', help='Path to the gtf folder', required=True)
p_ms.add_argument('--output', metavar='<String>', help='Path to the output folder', default="./mosyn_result/")
p_ms.add_argument('--window', metavar='<String>', help='The window size', default=0.1, type=float)
p_ms.add_argument("--complete", default=False, action='store_true', help="Complete synteny")
p_ms.add_argument('--mid', metavar='<String>', help='Binding sites ID', default="BS")
p_ms.add_argument('--idx', metavar='<String>', help='The starting number of the ID index', default=0, type=int)
p_ms.add_argument('--name', metavar='<String>', help='The Motif name', default="MOTIF")
p_ms.set_defaults(func=mosyn)

p_ra = subp.add_parser('runall', help='Run full pipeline')

p_ra.add_argument('--genome', metavar='<String>', help='Path to the genome folder', required=True)
p_ra.add_argument('--proteome', metavar='<String>', help='Path to the proteome folder', required=True)
p_ra.add_argument('--gtf', metavar='<String>', help='Path to the gtf folder', required=True)
p_ra.add_argument("--idconv", metavar="<String>", help="ID conversion folder", required=True)
p_ra.add_argument("--param", metavar="<String>", help="i-ADHoRe parameter file", required=True)
p_ra.add_argument("--pwm", metavar="<String>", help="PWM folder", required=True)

p_ra.add_argument('--output', metavar='<String>', help='Path to the output folder', default="./pipeline_result/")

p_ra.add_argument("--python2", metavar='<String>', help="Python 2 environment", required=True)
p_ra.add_argument("--ofopt", default=None, metavar="<String>",
                  help="Options for OrthoFinder. "
                       "Must be put inside a single quote, for example '-S diamond'. "
                       "Please refer to the OrthoFinder Manual.")

p_ra.add_argument("--pcol", default=0, metavar="<Integer>", help="Protein column")
p_ra.add_argument("--gcol", default=1, metavar="<Integer>", help="Gene column")
p_ra.add_argument("--csep", default="\t", metavar="<String>", help="Column separator")

p_ra.add_argument("--stopt", default=None, metavar="<String>",
                  help="Options for CREAD STORM. "
                       "Must be put inside a single quote, for example '-t 17'. "
                       "Please refer to the CREAD STORM Manual.")
p_ra.add_argument("--bcomp", default=False, action='store_true', help="Calculate base composition")

p_ra.add_argument('--window', metavar='<String>', help='The window size', default=0.1, type=float)
p_ra.add_argument("--complete", default=False, action='store_true', help="Complete synteny")
p_ra.add_argument('--mid', metavar='<String>', help='Binding sites ID', default="BS")
p_ra.add_argument('--idx', metavar='<String>', help='The starting number of the ID index', default=0, type=int)
p_ra.add_argument('--name', metavar='<String>', help='The Motif name', default="MOTIF")
p_ra.set_defaults(func=runall)


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
    elif argv[1] == 'mosyn':
        p_ms.print_help()
    elif argv[1] == 'runall':
        p_ra.print_help()
    exit(0)

args = p.parse_args()

args.func(args)
