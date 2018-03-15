"""
Generating summary
"""

import glob
import yaml


from misc.string import check_folder_path


def generate_orthofinder_summary(infolder, genus_index=-3, alignment_index=-2):
    """
    Generate OrthoFinder summary
    :param infolder: Input folder containing result
    :param genus_index: Genus index folder relative to result file
    :param alignment_index: Alignment index folder relative to result file
    :return:
    """

    outfile = "orthofinder_summary.csv"

    fout = open(outfile, 'w')

    print("Genus", "Alignment", "Number_of_Genes", "Number_of_Genes_in_Orthogroups", "Number_of_Unassigned_Genes",
          "Number_of_Orthogroups", sep=",", file=fout)

    infolder = check_folder_path(infolder)
    for summary in glob.glob(infolder+'**/Statistics_Overall.csv', recursive=True):

        path_elem = summary.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]

        fin = open(summary, 'r')
        lines = fin.readlines()
        fin.close()

        summary_keywords = ["Number of genes", "Number of genes in orthogroups", "Number of unassigned genes",
                            "Number of orthogroups"]
        summary_values = []

        for su in summary_keywords:
            for line in lines:
                line_elem = line.strip().split('\t')
                if line_elem[0] == su:
                    summary_values.append(int(line_elem[-1]))

        summary_write = [str(s) for s in summary_values]
        print(genus, alignment, ",".join(summary_write), sep=",", file=fout)

    fout.close()


def generate_mosyn_summary(infolder, genus_index=-5, alignment_index=-4, score_index=-3, pwm_index=-2):
    """
    Generate MoSyn summary
    :param infolder: Input folder containing result
    :param genus_index: Genus index folder relative to result file
    :param alignment_index: Alignment index folder relative to result file
    :param score_index: Score index folder relative to result file
    :param pwm_index: PWM index folder relative to result file
    :return:
    """

    outfile = "mosyn_detail_summary.csv"
    outfile0 = "mosyn_short_summary.csv"

    fout = open(outfile, 'w')
    fout0 = open(outfile0, 'w')

    print("Genus", "Alignment", "Score", "PWM", "Synteny_ID", "Number_of_Segments", "Number_of_Genes", "Number_of_CTCF",
          "Total_Length", sep=",", file=fout)

    print("Genus", "Alignment", "Score", "PWM", "Number_of_Synteny", "Number_of_Genes_in_Synteny",
          "Number_of_Synteny_containing_Motifs", "Number_of_Motifs_in_Synteny",
          "Average_Number_of_Genes_per_Synteny_per_Species", "Average_Length_per_Synteny_per_Species",
          "Average_Number_of_Motifs_per_Synteny_per_Species", sep=",", file=fout0)

    infolder = check_folder_path(infolder)
    for synt in glob.glob(infolder + '**/synteny.yaml', recursive=True):

        path_elem = synt.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]
        score = path_elem[score_index]
        pwm = path_elem[pwm_index]

        with open(synt, 'r') as stream:
            iadhore_dict_with_motifs = yaml.load(stream)

        synteny_length = 0
        synteny_genes = 0
        num_of_mult = 0

        synteny_motifs = 0
        synteny_contain = 0

        set_of_genes = set()
        set_of_motifs = set()

        for key, value in sorted(iadhore_dict_with_motifs.items()):

            num_of_mult += 1

            position_keys = sorted([k for k in value.keys() if k != "loops"])
            avg_mult_genes = len(position_keys)

            segment_keys = sorted([k for k in value[position_keys[0]].keys() if k != "motifs"])
            num_of_segments = len(segment_keys)

            mult_length = 0
            for sk in segment_keys:
                s_loc = [value[position_keys[0]][sk]["start"], value[position_keys[0]][sk]["end"],
                         value[position_keys[-1]][sk]["start"], value[position_keys[-1]][sk]["end"]]
                s_start = min(s_loc)
                s_end = max(s_loc)
                s_length = abs(s_end - s_start)
                mult_length += s_length
            avg_mult_length = mult_length / num_of_segments

            mult_motifs = 0
            mult_genes = 0

            avg_mult_motifs = 0

            check_contain = False

            for ke, val in value.items():

                if ke == "loops":
                    continue

                for k, v in val.items():

                    if k == "motifs":

                        if not check_contain:
                            check_contain = True

                        avg_mult_motifs += len(v)
                        for pair in v:
                            for motif in pair:
                                mult_motifs += 1
                                set_of_motifs.add(motif["motif"])

                    else:
                        mult_genes += 1
                        set_of_genes.add(v["gene"])

            alternative_mult_genes = avg_mult_genes * num_of_segments
            alternative_mult_motifs = avg_mult_motifs * num_of_segments

            assert alternative_mult_genes == mult_genes
            assert alternative_mult_motifs == mult_motifs

            if check_contain:
                synteny_contain += 1

            synteny_length += avg_mult_length
            synteny_genes += avg_mult_genes
            synteny_motifs += avg_mult_motifs

            print(genus, alignment, score, pwm, key, num_of_segments, mult_genes, mult_motifs, mult_length,
                  sep=",", file=fout)

        nr_genes = len(set_of_genes)
        nr_motifs = len(set_of_motifs)

        avg_synteny_length = synteny_length / num_of_mult
        avg_synteny_genes = synteny_genes / num_of_mult
        avg_synteny_motifs = synteny_motifs / synteny_contain

        print(genus, alignment, score, pwm, num_of_mult, nr_genes, synteny_contain, nr_motifs,
              avg_synteny_genes, avg_synteny_length, avg_synteny_motifs, sep=",", file=fout0)

    fout.close()
    fout0.close()
