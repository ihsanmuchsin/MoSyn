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


def generate_loop_summary(infolder, outfolder, genus_index=-5, alignment_index=-4, score_index=-3, pwm_index=-2):
    """
    Generate Loop summary
    :param infolder: Input folder containing result
    :param outfolder: Output folder
    :param genus_index: Genus index folder relative to result file
    :param alignment_index: Alignment index folder relative to result file
    :param score_index: Score index folder relative to result file
    :param pwm_index: PWM index folder relative to result file
    :return:
    """

    outfolder = check_folder_path(outfolder, True)

    outfile = outfolder + "loop_detail_summary.csv"
    outfile0 = outfolder + "loop_short_summary.csv"

    fout = open(outfile, 'w')
    fout0 = open(outfile0, 'w')

    print("Genus", "Alignment", "Score", "PWM", "Loop_ID", "Number_of_Segments", "Number_of_Genes", "Number_of_CTCF",
          "Total_Length", sep=",", file=fout)

    print("Genus", "Alignment", "Score", "PWM", "Number_of_Loops",
          "Number_of_Genes_in_Loops", "Number_of_Motifs_in_Loops",
          "Average_Number_of_Genes_per_Loops_per_Species", "Average_Length_per_Synteny_per_Species",
          "Average_Number_of_Motifs_per_Loops_per_Species", sep=",", file=fout0)

    infolder = check_folder_path(infolder)
    for loop_file in glob.glob(infolder + '**/loops.yaml', recursive=True):

        path_elem = loop_file.split('/')
        genus = path_elem[genus_index]
        alignment = path_elem[alignment_index]
        score = path_elem[score_index]
        pwm = path_elem[pwm_index]

        with open(loop_file, 'r') as stream:
            iadhore_dict_with_loops = yaml.load(stream)

        check_loops = False
        for key, value in sorted(iadhore_dict_with_loops.items()):
            if "loops" in value.keys():
                check_loops = True
                break

        if not check_loops:
            continue

        this_outdir = outfolder + "/".join([genus, alignment, score, pwm])
        this_outdir = check_folder_path(this_outdir, True)

        genes_and_motifs = this_outdir + "genes_and_motifs.txt"
        genes_only = this_outdir + "genes.txt"
        motifs_only = this_outdir + "motifs.txt"

        f_all = open(genes_and_motifs, 'w')
        f_gene = open(genes_only, 'w')
        f_mot = open(motifs_only, 'w')

        overall_length = 0
        overall_genes = 0
        num_of_loops = 0

        overall_motifs = 0

        set_of_genes = set()
        set_of_motifs = set()

        loop_index = 1
        for key, value in sorted(iadhore_dict_with_loops.items()):

            if "loops" not in value.keys():
                continue

            position_keys = sorted([k for k in value.keys() if k != "loops"])
            for loop in value["loops"]:

                num_of_loops += 1

                first_motif = loop["first"]
                last_motif = loop["last"]

                first_pos = loop["first_pos"]
                last_pos = loop["last_pos"]

                gc_keys = sorted([(m["genome"], m["chromosome"]) for m in first_motif])
                num_of_segments = len(gc_keys)

                genome_string = ",".join([str(g[0]) for g in gc_keys])
                chromosome_string = ",".join([str(g[-1]) for g in gc_keys])

                print("#loop_id=" + str(loop_index) + ";", "genome=" + genome_string + ";",
                      "chromosome=" + chromosome_string + ";", file=f_all)
                print("#loop_id=" + str(loop_index) + ";", "genome=" + genome_string + ";",
                      "chromosome=" + chromosome_string + ";", file=f_gene)
                print("#loop_id=" + str(loop_index) + ";", "genome=" + genome_string + ";",
                      "chromosome=" + chromosome_string + ";", file=f_mot)

                loop_loc = [(f["start"], f["end"], l["start"], l["end"]) for f, l in zip(first_motif, last_motif)]
                loop_length = sum([abs(max(l) - min(l)) for l in loop_loc])
                avg_loop_length = loop_length / num_of_segments

                check_start = False

                loop_genes = 0
                loop_motifs = 0

                for pk in position_keys:

                    if first_pos <= pk <= last_pos:

                        val = value[pk]
                        segment_keys = [s for s in val.keys() if s != "motifs"]

                        this_position_genes = []
                        this_gcount = 0
                        this_sgenes = []
                        for gk in gc_keys:
                            this_gene = None
                            for sk in segment_keys:
                                this_gk = (val[sk]["genome"], val[sk]["chromosome"])
                                if this_gk == gk:
                                    this_gene = val[sk]
                            if this_gene:
                                this_gcount += 1
                                this_sgenes.append(this_gene["gene"])
                                this_position_genes.append(this_gene["gene"] + this_gene["strand"])
                            else:
                                this_position_genes.append("-")

                        if this_position_genes and check_start:
                            loop_genes += this_gcount
                            set_of_genes.update(this_sgenes)
                            print("\t".join(this_position_genes), file=f_all)
                            print("\t".join(this_position_genes), file=f_gene)

                        if "motifs" not in val.keys():
                            continue

                        for pair in val["motifs"]:

                            if pair == first_motif:
                                check_start = True

                            this_position_motifs = []
                            this_mcount = 0
                            this_smot = []
                            for gk in gc_keys:
                                this_motif = None
                                for m in pair:
                                    m_gk = (m["genome"], m["chromosome"])

                                    if m_gk == gk:
                                        this_motif = m

                                if this_motif:
                                    this_mcount += 1
                                    this_smot.append(this_motif["motif"])
                                    this_position_motifs.append(this_motif["motif"] + this_motif["strand"])
                                else:
                                    this_position_motifs.append("-")

                            if this_position_motifs and check_start:
                                loop_motifs += this_mcount
                                set_of_motifs.update(this_smot)
                                print("\t".join(this_position_motifs), file=f_all)
                                print("\t".join(this_position_motifs), file=f_mot)

                            if pair == last_motif:
                                check_start = False

                avg_loop_genes = loop_genes / num_of_segments
                avg_loop_motifs = loop_motifs / num_of_segments

                overall_length += avg_loop_length
                overall_genes += avg_loop_genes
                overall_motifs += avg_loop_motifs

                print(genus, alignment, score, pwm, key, num_of_segments, loop_genes, loop_motifs, loop_length,
                      sep=",", file=fout)

                loop_index += 1

        f_all.close()
        f_gene.close()
        f_mot.close()

        nr_genes = len(set_of_genes)
        nr_motifs = len(set_of_motifs)

        avg_synteny_length = overall_length / num_of_loops
        avg_synteny_genes = overall_genes / num_of_loops
        avg_synteny_motifs = overall_motifs / num_of_loops

        print(genus, alignment, score, pwm, num_of_loops, nr_genes, nr_motifs,
              avg_synteny_genes, avg_synteny_length, avg_synteny_motifs, sep=",", file=fout0)

    fout.close()
    fout0.close()