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
