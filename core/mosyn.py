"""
MoSyn function
"""

import prep.gtf as pg


def restructure_gtf_dict(gtf_dict):
    """
    Restructure GTF dict. Only takes genes. Genome -> Chromosome -> Gene ID
    :param gtf_dict: Normal GTF dict
    :return: Restructured GTF dict
    """

    restructured_dict = dict()
    for key, value in gtf_dict.items():
        genome_dict = dict()
        for k, v in value.items():
            for v0 in v:
                if v0["feature"] == "gene":
                    gene_id = v0["attribute"]["gene_id"]
                    if k not in genome_dict.keys():
                        genome_dict[k] = dict()
                    genome_dict[k][gene_id] = v0

        restructured_dict[key] = genome_dict
    return restructured_dict


def add_location_to_iadhore_synteny(iadhore_synteny_dict, gtf_folder):
    """
    Add location info to iadhore synteny dictionary
    :param iadhore_synteny_dict: iadhore synteny dictionary
    :param gtf_folder: GTF folder
    :return: i-ADHoRe dictionary with location
    """

    gtf_dict = pg.gtf_to_dict_folder(gtf_folder)
    gtf_dict = restructure_gtf_dict(gtf_dict)

    for key, value in iadhore_synteny_dict.items():
        for k, v in value["segments"].items():

            genome = v["genome"]
            gene_start = v["first"]
            gene_end = v["last"]
            chromosome = str(v["list"])

            v["start"] = int(gtf_dict[genome][chromosome][gene_start]["start"])
            v["end"] = int(gtf_dict[genome][chromosome][gene_end]["end"])

            for k0, v0 in v["elements"].items():

                gene_id = v0["gene"]

                v0["start"] = int(gtf_dict[genome][chromosome][gene_id]["start"])
                v0["end"] = int(gtf_dict[genome][chromosome][gene_id]["end"])
                v0["strand"] = gtf_dict[genome][chromosome][gene_id]["strand"]

            first_element_location = list(v["elements"].values())[0]["start"]
            last_element_location = list(v["elements"].values())[-1]["start"]

            if first_element_location > last_element_location:
                v["inverse"] = True
            else:
                v["inverse"] = False

    return iadhore_synteny_dict


def add_motifs_into_synteny(iadhore_synteny_dict, motifs_gtf_folder):
    """
    Add motifs as a list of set(motif_id, distance_to_element) as a value of element["motifs"] in iadhore synteny dict
    :param iadhore_synteny_dict: i-ADHoRe synteny dict
    :param motifs_gtf_folder: motifs in GTF format
    :return: i-ADHoRe synteny with added motifs
    """

    motifs_dict = pg.gtf_to_dict_folder(motifs_gtf_folder)

    for key, value in motifs_dict.items():
        for k, v in value.items():
            for v1 in v:

                for key0, value0 in iadhore_synteny_dict.items():
                    for k0, v0 in value0["segments"].items():

                        genome = v0["genome"]
                        chromosome = str(v0["list"])

                        if not (genome == key and chromosome == k):
                            continue

                        motif_start = int(v1["start"])
                        motif_end = int(v1["end"])
                        motif_strand = v1["strand"]
                        segment_start = v0["start"]
                        segment_end = v0["end"]

                        if not (motif_start >= segment_start and motif_end <= segment_end):
                            continue

                        inverse = v0["inverse"]
                        element_keys = list(v0["elements"].keys())

                        for i in range(len(element_keys)-1):

                            k20 = element_keys[i]
                            k21 = element_keys[i+1]

                            this_element = v0["elements"][k20]
                            this_element_start = this_element["start"]
                            this_element_end = this_element["end"]

                            next_element = v0["elements"][k21]
                            next_element_start = next_element["start"]
                            next_element_end = next_element["end"]

                            check_motif = False
                            if inverse:
                                if motif_end <= this_element_end and motif_start > next_element_end:
                                    check_motif = True
                                    distance = abs(this_element_end - motif_end)
                            else:
                                if motif_start >= this_element_start and motif_end < next_element_start:
                                    check_motif = True
                                    distance = abs(motif_start - this_element_start)

                            if check_motif:
                                if "motifs" not in this_element.keys():
                                    this_element["motifs"] = []
                                motif_id = v1["attribute"]["motif_id"]
                                this_dict = {"motif": motif_id,
                                             "start": motif_start,
                                             "end": motif_end,
                                             "strand": motif_strand,
                                             "distance": distance,
                                             "chromosome": chromosome,
                                             "genome": genome}
                                this_element["motifs"].append(this_dict)

    return iadhore_synteny_dict


def restructure_iadhore_dict_to_position(iadhore_synteny_dict):
    """
    Change the structure the i-ADHoRe dictionary to Multiplicon -> Position -> Genome -> Element.
    The resulted dictionary only retain ID, chromosome, start, end, strand, and motifs
    :param iadhore_synteny_dict: Normal i-ADHoRe synteny dict
    :return: Restructured i-ADHoRe synteny dict
    """

    restructured_dict = dict()
    for key, value in iadhore_synteny_dict.items():
        position = dict()
        for k, v in value["segments"].items():
            genome = v["genome"]
            chromosome = v["list"]
            for k0, v0 in v["elements"].items():

                pos = int(v0["position"])
                if pos not in position.keys():
                    position[pos] = dict()

                gene_dict = dict()
                gene_dict["start"] = v0["start"]
                gene_dict["end"] = v0["end"]
                gene_dict["strand"] = v0["strand"]
                gene_dict["gene"] = v0["gene"]
                gene_dict["chromosome"] = chromosome
                gene_dict["genome"] = genome

                try:
                    gene_dict["motifs"] = v0["motifs"]
                except KeyError:
                    pass

                position[pos][k] = gene_dict

        restructured_dict[key] = position

    return restructured_dict


def align_motifs_in_synteny(iadhore_position_with_motifs_dict, window=0.1):
    """
    Align motifs into i-ADHoRe synteny
    :param iadhore_position_with_motifs_dict: i-ADHoRe position with motifs
    :param window: Window size (percentage)
    :return: i-ADHoRe position with aligned motifs
    """

    for key, value in iadhore_position_with_motifs_dict.items():
        for ke, val in value.items():

            all_pairs = []

            for k in val.keys():

                if "motifs" not in val[k].keys():
                    continue

                other_keys = set(val.keys())
                other_keys.discard(k)

                for motif in val[k]["motifs"]:

                    dist = motif["distance"]
                    this_strand = val[k]["strand"] == motif["strand"]
                    pairs = []

                    for k0 in other_keys:

                        if "motifs" not in val[k0].keys():
                            continue

                        selected_pair = []

                        for mot in val[k0]["motifs"]:

                            that_strand = val[k0]["strand"] == mot["strand"]

                            if not(this_strand == that_strand):
                                continue

                            d = mot["distance"]
                            diff = abs(dist - d)

                            if diff <= (window*dist):
                                selected_pair.append(mot)
                                continue

                        if not selected_pair:
                            continue

                        if not pairs:
                            for sl in selected_pair:
                                pairs.append([sl])
                            continue

                        mod_pairs = []
                        for pair in pairs:
                            for sl in selected_pair:
                                mod_pairs.append(pair + [sl])
                                all_pairs.append(pair + [sl])
                        pairs = mod_pairs

                    if not pairs:
                        continue

                    for pair in pairs:
                        all_pairs.append(pair + [motif])

            # modified all pairs
            nr_all_pairs = []
            num_of_species = len(val)
            for i in range(num_of_species, 1, -1):
                kmer_pairs = dict()
                for el in all_pairs:
                    if len(el) == i:
                        el_dist = [m["distance"] for m in el]
                        el_sum = sum(list(el_dist))

                        if el_sum not in kmer_pairs.keys():
                            kmer_pairs[el_sum] = []
                        kmer_pairs[el_sum].append(el)

                for su, ls in sorted(kmer_pairs.items()):
                    for l in ls:
                        if not nr_all_pairs:
                            nr_all_pairs.append(l)
                            continue

                        check_nr = True
                        for nr in nr_all_pairs:
                            for e in l:
                                if e in nr:
                                    check_nr = False
                                    break
                            if not check_nr:
                                break
                        if check_nr:
                            nr_all_pairs.append(l)

            # then add singleton
            single = []
            for k in val.keys():

                if "motifs" not in val[k].keys():
                    continue

                for motif in val[k]["motifs"]:
                    check_single = True
                    for pair in nr_all_pairs:
                        if motif in pair:
                            check_single = False
                            break
                    if check_single:
                        single.append([motif])

            nr_final_pairs = nr_all_pairs + single

            # delete working motifs
            for k, v in val.items():

                if "motifs" not in v.keys():
                    continue

                del v["motifs"]

            # append final result
            if nr_final_pairs:
                val["motifs"] = nr_final_pairs

    return iadhore_position_with_motifs_dict


def get_complete_motifs_synteny(iadhore_with_aligned_motifs_dict):
    """
    Get complete pairs, which have pair in all of the segments
    :param iadhore_with_aligned_motifs_dict: i-ADHoRe dict with aligned motifs
    :return: complete motifs synteny
    """

    for key, value in iadhore_with_aligned_motifs_dict.items():
        for ke, val in value.items():

            if "motifs" not in val.keys():
                continue

            num_of_species = len(list(val.keys())) - 1

            complete_pairs = []
            for pairs in val["motifs"]:

                if len(pairs) < num_of_species:
                    continue

                complete_pairs.append(pairs)

            if not complete_pairs:
                del val["motifs"]
                continue

            val["motifs"] = complete_pairs

    return iadhore_with_aligned_motifs_dict
