import os

from utils.classification.pgdb_filter_test import process_pgdb_folder

savi_type_list = ['AIPP', 'CAPP', 'CVP', 'MANUAL', 'NPP', 'UPP', 'problem']


def read_rxn_list(rxn_list_path):
    rxn_list = set()
    with open(rxn_list_path, 'r') as op:
        for line in op:
            line = line.strip()
            if line != '':
                rxn_list.add(line)
    return rxn_list


def subset_dict(input_dict, key_list):
    output_dict = {}
    for key in sorted(key_list):
        try:
            output_dict.setdefault(key, input_dict[key])
        except KeyError:
            continue
    return output_dict


def get_vals_of_dict(input_dict):
    dict_vals = set()
    for key in input_dict:
        vals = input_dict[key]
        if type(vals) is list or type(vals) is set:
            dict_vals.update(vals)
        else:
            dict_vals.add(vals)
    return dict_vals


def reverse_dict(input_dict):
    output_dict = {}
    for key in input_dict:
        vals = input_dict[key]
        if type(vals) is list or type(vals) is set:
            for val in sorted(vals):
                try:
                    output_dict[val].add(key)
                except KeyError:
                    output_dict.setdefault(val, {key})
        else:
            try:
                output_dict[vals].add(key)
            except KeyError:
                output_dict.setdefault(vals, {key})
    return output_dict


def read_savi_types(savi_type_path):
    pathway_to_savi = {}
    with open(savi_type_path, 'r') as stp:
        for line in stp:
            info = [i.strip() for i in line.split('\t')]
            try:
                pathway_to_savi[info[0]].add(info[1])
            except (KeyError, IndexError) as e:
                pathway_to_savi.setdefault(info[0], {info[1]})
    return pathway_to_savi


def partition_reaction_to_pathways(op, prefix, reaction_to_pathway_before, reaction_to_pathway_after, pathway_to_savi,
                                   suffix='\n'):
    reaction_before = sorted(reaction_to_pathway_before.keys())
    reaction_after = sorted(reaction_to_pathway_after.keys())

    pathway_before = sorted(get_vals_of_dict(reaction_to_pathway_before))
    pathway_after = sorted(get_vals_of_dict(reaction_to_pathway_after))

    pathway_types_before = reverse_dict(subset_dict(pathway_to_savi, pathway_before))
    pathway_types_after = reverse_dict(subset_dict(pathway_to_savi, pathway_after))

    count_reaction_before = len(reaction_before)
    count_reaction_after = len(reaction_after)
    precent_reaction_after = float(count_reaction_after) / float(count_reaction_before)
    count_reaction_removed = count_reaction_before - count_reaction_after
    percent_reaction_removed = float(count_reaction_removed) / float(count_reaction_before)
    count_pathway_before = len(pathway_before)
    count_pathway_after = len(pathway_after)
    precent_pathway_after = float(count_pathway_after) / float(count_pathway_before)
    count_pathway_removed = count_pathway_before - count_pathway_after
    percent_pathway_removed = float(count_pathway_removed) / float(count_pathway_before)

    output_list = [count_reaction_before, count_reaction_after, precent_reaction_after, count_reaction_removed,
                   percent_reaction_removed, count_pathway_before, count_pathway_after, precent_pathway_after,
                   count_pathway_removed, percent_pathway_removed]
    for t in sorted(savi_type_list):
        try:
            count_type_before = len(pathway_types_before[t])
        except KeyError:
            count_type_before = 0
        try:
            count_type_after = len(pathway_types_after[t])
        except KeyError:
            count_type_after = 0
        count_type_removed = count_type_before - count_type_after
        try:
            percent_type_after = float(count_type_after) / float(count_type_before)
        except ZeroDivisionError:
            percent_type_after = 0
        try:
            percent_type_removed = float(count_type_removed) / float(count_type_before)
        except ZeroDivisionError:
            percent_type_removed = 0
        output_list += [count_type_before, count_type_after, percent_type_after, count_type_removed, percent_type_removed]

    op.write(prefix + '\t' + '\t'.join([str(o) for o in output_list]) + suffix)
    return reaction_before, reaction_after, pathway_before, pathway_after, pathway_types_before, pathway_types_after


def process_pgdbs(pre_savi_folder, release_folder, rxn_list, pathway_to_savi, output_path):
    with open(output_path, 'w') as op:
        op.write('\t'.join(['PGDB', 'Type', 'Reaction', 'ReactionKept', 'ReactionKeptPercentage',
                            'ReactionRemoved', 'ReactionRemovedPercentage',
                            'Pathway', 'PathwayKept', 'PathwayKeptPercentage', 'PathwayRemoved',
                            'PathwayRemovedPercentage',
                            'PathwayAIPP', 'PathwayAIPPKept', 'PathwayAIPPKeptPercentage', 'PathwayAIPPRemoved', 'PathwayAIPPRemovedPercentage',
                            'PathwayCAPP', 'PathwayCAPPKept', 'PathwayCAPPKeptPercentage', 'PathwayCAPPRemoved', 'PathwayCAPPRemovedPercentage',
                            'PathwayCVP', 'PathwayCVPKept', 'PathwayCVPKeptPercentage', 'PathwayCVPRemoved', 'PathwayCVPRemovedPercentage',
                            'PathwayMANUAL', 'PathwayMANUALKept', 'PathwayMANUALKeptPercentage', 'PathwayMANUALRemoved', 'PathwayMANUALRemovedPercentage',
                            'PathwayNPP', 'PathwayNPPKept', 'PathwayNPPKeptPercentage', 'PathwayNPPRemoved', 'PathwayNPPRemovedPercentage',
                            'PathwayUPP', 'PathwayUPPKept', 'PathwayUPPKeptPercentage', 'PathwayUPPRemoved', 'PathwayUPPRemovedPercentage',
                            'Pathwayproblem', 'PathwayproblemKept', 'PathwayproblemKeptPercentage', 'PathwayproblemRemoved', 'PathwayproblemRemovedPercentage']) + '\n')
        for f in sorted(set(os.listdir(pre_savi_folder)) & set(os.listdir(release_folder))):
            if f.startswith('.'):
                continue
            pre_savi_pgdb = os.path.join(pre_savi_folder, f)
            _, pre_savi_pathways_dat_reaction_to_pathway, _ = process_pgdb_folder(pre_savi_pgdb)
            pre_savi_pathways_dat_reaction_to_pathway_subset = subset_dict(pre_savi_pathways_dat_reaction_to_pathway,
                                                                           rxn_list)
            partition_reaction_to_pathways(op, f + '\tPreSAVI', pre_savi_pathways_dat_reaction_to_pathway,
                                           pre_savi_pathways_dat_reaction_to_pathway_subset, pathway_to_savi)

            release_pgdb = os.path.join(release_folder, f)
            _, release_pathways_dat_reaction_to_pathway, _ = process_pgdb_folder(release_pgdb)
            release_pathways_dat_reaction_to_pathway_subset = subset_dict(release_pathways_dat_reaction_to_pathway,
                                                                          rxn_list)
            partition_reaction_to_pathways(op, f + '\tRelease', release_pathways_dat_reaction_to_pathway,
                                           release_pathways_dat_reaction_to_pathway_subset, pathway_to_savi)


if __name__ == '__main__':
    fscore_rxn_list = read_rxn_list('/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/rxn.fscore_0.5.txt')
    precision_rxn_list = read_rxn_list(
        '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/rxn.precision_0.5.txt')
    size_rxn_list = read_rxn_list('/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/rxn.size_10.txt')

    filter_test_output_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/filtering_wo_savi'
    pre_savi_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_preSAVI_02_16_21'
    release_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_release_02_16_21'

    savi_types = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/savi.types.txt'
    pathway_to_savi = read_savi_types(savi_types)
    # print(pathway_to_savi)
    process_pgdbs(pre_savi_folder, release_folder, fscore_rxn_list, pathway_to_savi,
                  os.path.join(filter_test_output_folder, 'fscore.tsv'))

    process_pgdbs(pre_savi_folder, release_folder, precision_rxn_list, pathway_to_savi,
                  os.path.join(filter_test_output_folder, 'precision.tsv'))

    process_pgdbs(pre_savi_folder, release_folder, size_rxn_list, pathway_to_savi,
                  os.path.join(filter_test_output_folder, 'size.tsv'))
