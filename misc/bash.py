"""
Create BASH for automation
"""

from misc.string import check_folder_path


def create_storm_bash(script_file, genome_dir, pwm_dir, outdir, bash_outfile,
                      min_threshold, max_threshold, increment):

    outdir = check_folder_path(outdir)

    fout = open(bash_outfile, 'w')
    for i in range(min_threshold, max_threshold+1, increment):

        storm_opt = '\'-t ' + str(i) + '\''
        sub_outdir = outdir + str(i)

        print('python', script_file, 'storm', '--input', genome_dir, '--pwm', pwm_dir, '--output', sub_outdir,
              '--opt', storm_opt, '--bcomp', file=fout)
    fout.close()


def create_mosyn_bash(result_folder, material_folder, output_folder, pwm_names,
                      list_genus, list_alignment, range_score, script_file, bash_outfile):

    result_folder = check_folder_path(result_folder)
    material_folder = check_folder_path(material_folder)
    output_folder = check_folder_path(output_folder)

    fout = open(bash_outfile, 'w')

    for genus in list_genus:
        for alignment in list_alignment:
            for score in range_score:
                for pwm in pwm_names:

                    iadhore_output = result_folder + "/".join(["IADHORE", genus, alignment])
                    iadhore_output = check_folder_path(iadhore_output)

                    storm_output = result_folder + "/".join(["STORM", genus, str(score), pwm])
                    storm_output = check_folder_path(storm_output)

                    gtf_folder = material_folder + "/".join([genus, "GTF"])
                    gtf_folder = check_folder_path(gtf_folder)

                    this_outdir = output_folder + "/".join([genus, alignment, str(score), pwm])
                    this_outdir = check_folder_path(this_outdir)

                    print("python", script_file, "mosyn", "--iadhore", iadhore_output, "--storm", storm_output,
                          "--gtf", gtf_folder, "--output", this_outdir,
                          "--complete", "--mid", "CTCF", "--idx", 0, "--name", "CTCF", file=fout)

    fout.close()


def create_loop_bash(result_folder, output_folder, pwm_names, list_genus, list_alignment, range_score,
                     script_file, bash_outfile, min_length=80000, max_length=800000):

    result_folder = check_folder_path(result_folder)
    output_folder = check_folder_path(output_folder)

    fout = open(bash_outfile, 'w')

    for genus in list_genus:
        for alignment in list_alignment:
            for score in range_score:
                for pwm in pwm_names:

                    this_input = result_folder + "/".join([genus, alignment, str(score), pwm])
                    this_input = check_folder_path(this_input)
                    this_input += "synteny.yaml"

                    this_outdir = output_folder + "/".join([genus, alignment, str(score), pwm])
                    this_outdir = check_folder_path(this_outdir)

                    print("python", script_file, "idloop", "--input", this_input, "--output", this_outdir,
                          "--min", min_length, "--max", max_length, file=fout)

    fout.close()



def generate_mosyn_short_summary(infolder, outfile, genus_index=-5, alignment_index=-4, score_index=-3, pwm_index=-2):
    """
    Generate detailed i-ADHoRe summary.
    Result structure Genus -> Alignment -> Result
    GTF structure Material -> Genus -> GTF
    :param material_folder: Folder containing GTF folder
    :param infolder: Folder containing i-ADHoRe result
    :param outfile: Output folder
    :param genus_index: Genus folder index from result
    :param alignment_index: Alignment folder index from result
    :return:
    """

    material_folder = check_folder_path(material_folder)

    fout = open(outfile, 'w')

    print("Genus", "Alignment", "Number_of_Synteny", "Number_of_Genes_in_Synteny",
          "Average_Number_of_Genes_per_Synteny", "Average_Length_per_Synteny",
          sep=",", file=fout)

    infolder = check_folder_path(infolder)
    for synt in glob.glob(infolder + '**/synteny.yaml', recursive=True):

        path_elem = mult.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]
        pwm = path_elem[pwm_index]
        score = path_elem[score_index]

        with open(synt, 'r') as stream:
            iadhore_dict_with_motifs = yaml.load(stream)

        synteny_length = 0
        synteny_genes = 0
        num_of_mult = 0

        set_of_genes = set()

        for key, value in iadhore_dict_with_motifs.items():

            number_of_segments = 0
            mult_genes = 0
            mult_length = 0

            for ke, val in value["segments"].items():
                number_of_segments += 1
                length = abs(val["start"] - val["end"])
                mult_length += length
                mult_genes += len(val["elements"])

                for v in val["elements"].values():
                    set_of_genes.add(v["gene"])

            mult_avg_len = mult_length / number_of_segments
            mult_avg_genes = mult_genes / number_of_segments

            synteny_length += mult_avg_len
            synteny_genes += mult_avg_genes
            num_of_mult += 1

        nr_genes = len(set_of_genes)
        synteny_avg_length = synteny_length / num_of_mult
        synteny_avg_genes = synteny_genes / num_of_mult

        print(genus, alignment, num_of_mult, nr_genes, synteny_avg_genes, synteny_avg_length, sep=",", file=fout)

    fout.close()
