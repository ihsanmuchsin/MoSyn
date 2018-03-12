"""
MoSyn function
"""

import numpy
import json
import yaml
import pickle

import prep.gtf as pg
import prep.storm as ps
import prep.iadhore as pi
from misc.string import check_folder_path


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


def add_motifs_into_synteny(iadhore_location, motifs_gtf_folder):
    """
    Add motifs as a list of set(motif_id, distance_to_element) as a value of element["motifs"] in iadhore synteny dict
    :param iadhore_location: i-ADHoRe synteny dict
    :param motifs_gtf_folder: motifs in GTF format
    :return: i-ADHoRe synteny with added motifs
    """

    motifs_dict = pg.gtf_to_dict_folder(motifs_gtf_folder)

    for key, value in motifs_dict.items():
        for k, v in value.items():
            for v1 in v:

                for key0, value0 in iadhore_location.items():
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

    return iadhore_location


def restructure_iadhore_dict_to_position(iadhore_location):
    """
    Change the structure the i-ADHoRe dictionary to Multiplicon -> Position -> Genome -> Element.
    The resulted dictionary only retain ID, chromosome, start, end, strand, and motifs
    :param iadhore_location: Normal i-ADHoRe synteny dict
    :return: Restructured i-ADHoRe synteny dict
    """

    restructured_dict = dict()
    for key, value in iadhore_location.items():
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

            # sort final pairs
            mean_distance = dict()
            for pair in nr_final_pairs:
                pr_dist = [p["distance"] for p in pair]
                pr_mean = numpy.mean(pr_dist)
                mean_distance[pr_mean] = pair

            sorted_pairs = []
            for di, pr in sorted(mean_distance.items()):
                sorted_pairs.append(pr)

            # delete working motifs
            for k, v in val.items():

                if "motifs" not in v.keys():
                    continue

                del v["motifs"]

            # append final result
            if sorted_pairs:
                val["motifs"] = sorted_pairs

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


def dump_aligned_motifs_to_flat_synteny(iadhore_with_aligned_motifs, outfile):
    """
    Dump aligned motifs to synteny flat file
    :param iadhore_with_aligned_motifs: i-ADHoRe with aligned motifs
    :param outfile: synteny flat file
    :return:
    """

    fout = open(outfile, 'w')

    for key, value in sorted(iadhore_with_aligned_motifs.items()):

        genome_chromosome_list = []
        for ke, val in value.items():

            segment_keys = list(val.keys())
            if "motifs" in val.keys():
                segment_keys.remove("motifs")

            for k in segment_keys:
                genome_chromosome = [val[k]["genome"], val[k]["chromosome"]]
                if genome_chromosome not in genome_chromosome_list:
                    genome_chromosome_list.append(genome_chromosome)

        genome_chromosome_list = sorted(genome_chromosome_list)

        genome_str = ",".join([str(g[0]) for g in genome_chromosome_list])
        chromosome_str = ",".join([str(g[1]) for g in genome_chromosome_list])
        print("#synteny_id="+str(key)+";", "genome="+genome_str+";", "chromosome="+chromosome_str+";",
              file=fout)

        for ke, val in sorted(value.items()):
            this_position_genes = []
            for gn_chr in genome_chromosome_list:
                this_gene = None
                for k, v in val.items():

                    if k == "motifs":
                        continue

                    g_genome = v["genome"]
                    g_chromosome = v["chromosome"]

                    g_gc = [g_genome, g_chromosome]
                    if not (g_gc == gn_chr):
                        continue

                    this_gene = v
                    break

                if this_gene:
                    this_position_genes.append(this_gene["gene"]+this_gene["strand"])
                else:
                    this_position_genes.append("-")

            print("\t".join(this_position_genes), file=fout)

            if "motifs" not in val.keys():
                continue

            for pair in val["motifs"]:
                this_position_motifs = []
                for gc in genome_chromosome_list:
                    this_p = None
                    for p in pair:
                        p_genome = p["genome"]
                        p_chromosome = p["chromosome"]
                        p_gc = [p_genome, p_chromosome]

                        if not (p_gc == gc):
                            continue

                        this_p = p
                        break

                    if this_p:
                        this_position_motifs.append(this_p["motif"]+this_p["strand"])
                    else:
                        this_position_motifs.append("-")

                print("\t".join(this_position_motifs), file=fout)

    fout.close()


def dump_aligned_motifs_to_serial(iadhore_with_aligned_motifs, outfile, ftype='yaml'):
    """
    Dump aligned motifs to serialized file
    :param ftype: File type, e.g., YAML JSON, Pickle
    :param iadhore_with_aligned_motifs: i-ADHoRe with aligned motifs
    :param outfile: serialized file
    :return:
    """

    fout = open(outfile, 'w')

    if ftype == "json":
        json.dump(iadhore_with_aligned_motifs, fout, indent=5)
    elif ftype == "yaml":
        yaml.dump(iadhore_with_aligned_motifs, fout)
    elif ftype == "pickle":
        fout.close()
        fout = open(outfile, 'wb')
        pickle.dump(iadhore_with_aligned_motifs, fout)

    fout.close()


def run_mosyn(iadhore_output_folder, storm_output_folder, gtf_folder, outfolder,
              window=0.1, complete=True, binding_site_id="BS", id_start_index=0, motif_name="MOTIF"):
    """
    Run MoSyn
    :param motif_name: The name of the motif
    :param id_start_index: The binding site id start index
    :param binding_site_id: The binding site id
    :param window: Window size of alignment
    :param iadhore_output_folder: i-ADHoRe output folder
    :param storm_output_folder: CREAD STORM output folder
    :param gtf_folder: GTF folder
    :param outfolder: MoSyn result folder
    :param complete: Complete alignment
    :return:
    """

    outfolder = check_folder_path(outfolder, True)

    working_directory = outfolder + 'working_directory/'
    working_directory = check_folder_path(working_directory, True)

    motifs_gtf_folder = working_directory + 'motifs_gtf/'
    motifs_gtf_folder = check_folder_path(motifs_gtf_folder, True)
    ps.transfac_to_gtf_folder(storm_output_folder, motifs_gtf_folder, binding_site_id, id_start_index, motif_name)

    iadhore_dict = pi.iadhore_result_folder_to_dict(iadhore_output_folder)
    if complete:
        iadhore_dict = pi.get_complete_synteny_dict(iadhore_dict)
    iadhore_with_location = add_location_to_iadhore_synteny(iadhore_dict, gtf_folder)
    iadhore_with_motifs = add_motifs_into_synteny(iadhore_with_location, motifs_gtf_folder)
    restructured_iadhore = restructure_iadhore_dict_to_position(iadhore_with_motifs)
    iadhore_with_pairs = align_motifs_in_synteny(restructured_iadhore, window)
    if complete:
        iadhore_with_pairs = get_complete_motifs_synteny(iadhore_with_pairs)

    flat_synteny = outfolder + "synteny.txt"
    serial_synteny = outfolder + "synteny.yaml"

    dump_aligned_motifs_to_flat_synteny(iadhore_with_pairs, flat_synteny)
    dump_aligned_motifs_to_serial(iadhore_with_pairs, serial_synteny)
