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

