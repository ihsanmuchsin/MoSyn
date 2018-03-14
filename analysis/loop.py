"""
Loop identification
"""

import yaml
import numpy


def identify_loops_in_synteny(motifs_synteny_yaml_file, min_length=10000, max_length=1000000):
    """
    Identify k-mer to k-mer loops
    :param motifs_synteny_yaml_file: MoSyn YAML output
    :param min_length: The minimum length of the loops
    :param max_length: The maximum length of the loops
    :return: MoSyn dictionary with loops
    """

    with open(motifs_synteny_yaml_file, 'r') as stream:
        mosyn_dict = yaml.load(stream)

    for key, value in sorted(mosyn_dict.items()):

        key_counter = 0
        for ke, val in sorted(value.items()):

            if "motifs" in val.keys():
                key_counter += 1

        if key_counter < 2:
            continue

        loops_size = dict()
        for ke, val in sorted(value.items()):
            if "motifs" in val.keys():
                for pair in val["motifs"]:

                    current_first = pair
                    current_pos = ke
                    current_strand = {(m["genome"], m["chromosome"]): m["strand"] for m in pair}
                    current_last = None
                    current_size = 0
                    last_pos = None
                    current_gc = set([(m["genome"], m["chromosome"]) for m in pair])

                    for ke0, val0 in sorted(value.items()):
                        if "motifs" in val0.keys():
                            for pair0 in val0["motifs"]:

                                pos_distance = ke0 - current_pos

                                if pos_distance < 2:
                                    break

                                this_gc = set([(m["genome"], m["chromosome"]) for m in pair0])
                                if not (current_gc == this_gc):
                                    continue

                                this_strand = {(m["genome"], m["chromosome"]): m["strand"] for m in pair0}
                                if current_strand == this_strand:
                                    continue

                                start = [m["start"] for m in current_first]
                                end = [m["start"] for m in pair0]
                                size_list = [abs(x - y) for x, y in zip(start, end)]
                                avg_size = numpy.mean(size_list)

                                if min_length <= avg_size <= max_length:

                                    size = sum(size_list)

                                    if (not current_last) or (size > current_size):
                                        current_last = pair0
                                        current_size = size
                                        last_pos = ke0

                    if current_last:

                        loop = {
                            "first": current_first,
                            "last": current_last,
                            "first_pos": current_pos,
                            "last_pos": last_pos
                        }

                        loop_size = int(current_size)

                        if loop_size not in loops_size.keys():
                            loops_size[loop_size] = []

                        loops_size[loop_size].append(loop)

        if not loops_size:
            continue

        selected_loops = []
        for loop_size in sorted(loops_size.keys(), reverse=True):
            for loop in loops_size[loop_size]:

                if not selected_loops:
                    selected_loops.append(loop)
                    continue

                check_select = True
                for sl in selected_loops:

                    if loop["first_pos"] >= sl["first_pos"] and loop["last_pos"] <= sl["last_pos"]:
                        check_select = False
                        break

                if check_select:
                    selected_loops.append(loop)

        pos_loops = dict()
        for loop in selected_loops:
            pos_key = (loop["first_pos"], loop["last_pos"])
            pos_loops[pos_key] = loop

        loops = []
        for pos, loop in sorted(pos_loops.items()):
            loops.append(loop)

        if loops:
            mosyn_dict[key]["loops"] = loops

    return mosyn_dict


def print_loops_to_csv(mosyn_loops_dict, outfile):
    """
    Print MoSyn loops to csv file
    :param mosyn_loops_dict: MoSyn loops dictionary
    :param outfile: .csv file
    :return:
    """

    fout = open(outfile, 'w')

    start_index = 1
    for key, value in sorted(mosyn_loops_dict.items()):

        if "loops" not in value.keys():
            continue

        for loop in value["loops"]:

            print_value = []
            segments = sorted([(m["genome"], m["chromosome"]) for m in loop["first"]])
            for seg in segments:
                l_motifs = [m for m in loop["first"]+loop["last"] if (m["genome"], m["chromosome"]) == seg]
                l_start = [m["start"] for m in l_motifs]
                l_end = [m["end"] for m in l_motifs]
                l_loc = sorted(l_start+l_end)
                l_id = sorted([m["motif"] for m in l_motifs])
                print_value += [seg[0], seg[-1], l_id[0], l_id[-1], l_loc[0], l_loc[-1]]

            print_string = [str(p) for p in print_value]

            print(start_index, key, ",".join(print_string), sep=",", file=fout)
            start_index += 1

    fout.close()
