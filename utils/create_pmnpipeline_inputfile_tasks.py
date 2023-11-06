from argparse import ArgumentParser
import os
import sys

import itertools

import re


def read_dat_itr(fp, id_field, split_str):
    unique_id, attrs = None, []
    for line in fp:
        line = line.strip()
        l = line.split(split_str)
        if l[0] == id_field:
            if unique_id: yield (unique_id, ''.join(attrs))
            unique_id, attrs = line, []
        else:
            attrs.append(line + "\n")
    if unique_id: yield (unique_id, ''.join(attrs))


def get_versions_from_input(file_path):
    with open(file_path, 'r', encoding='iso-8859-1') as fp:
        for line in fp:
            if line.startswith('# Version:') or line.startswith('CC   Release of'):
                version = line.replace('# Version:', '').replace('CC   Release of', '').strip()
                return version


def get_ec_mapping_from_reactions(file_path, version, output_folder):
    version = str(version)
    # input_folder = os.path.dirname(file_path)
    rxn_ec_name = 'metacyc-RXN-EC.mapping'
    rxn_ec_path = os.path.join(output_folder, rxn_ec_name)
    rxn_official_ec_name = 'metacyc-RXN-official-EC.mapping'
    rxn_official_ec_path = os.path.join(output_folder, rxn_official_ec_name)
    sub_reaction_name = 'metacyc-sub-reactions'
    sub_reaction_path = os.path.join(output_folder, sub_reaction_name)

    with open(file_path, 'r', encoding='iso-8859-1') as ip:
        with open(rxn_ec_path, 'w', encoding='iso-8859-1') as rep, \
                open(rxn_official_ec_path, 'w', encoding='iso-8859-1') as roep, \
                open(sub_reaction_path, 'w', encoding='iso-8859-1') as srp:
            rep.write('# MetaCyc Version:\t' + version + '\n')
            roep.write('# MetaCyc Version:\t' + version + '\n')
            srp.write('# MetaCyc Version:\t' + version + '\n')
            for u_id, attrs in read_dat_itr(ip, "UNIQUE-ID", " - "):
                unique_id = u_id.replace('UNIQUE-ID -', '').strip()
                attributes = attrs.strip().split('\n')

                ec_numbers = set()
                reaction_list = set()
                official_ec_numbers = set()

                for index, entries in enumerate(attributes):
                    fields = entries.split(" - ")
                    if fields[0] == 'EC-NUMBER':
                        ec_num = fields[1].replace('EC-', '').strip()
                        ec_numbers.add(ec_num)
                        try:
                            check_if_official = attributes[index + 1].split(" - ")
                            if check_if_official[0] == '^OFFICIAL?' and check_if_official[1].strip() == 'T':
                                official_ec_numbers.add(ec_num)
                        except IndexError:
                            pass
                    if fields[0] == 'REACTION-LIST':
                        reaction_list.add(fields[1].strip())

                for ec_num in ec_numbers:
                    rep.write(unique_id + '\t' + ec_num + '\n')
                for official_ec_num in official_ec_numbers:
                    roep.write(official_ec_num + '\t' + unique_id + '\n')
                for rxn in reaction_list:
                    srp.write(unique_id + '\t' + rxn + '\n')


def get_metacyc_common_name(reaction_path, enrxn_path, version, output_folder):
    # input_folder = os.path.dirname(reaction_path)
    metacyc_common_name = "metacyc-rxn-name-mapping"
    metacyc_common_path = os.path.join(output_folder, metacyc_common_name)
    common_name_not_found = os.path.join(output_folder,
                                         metacyc_common_name + '.name_not_found.' + str(version) + '.txt')
    enrxn_common_name = {}
    with open(enrxn_path, 'r', encoding='iso-8859-1') as ep:
        for u_id, attrs in read_dat_itr(ep, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            attributes = attrs.strip().split('\n')

            for index, entries in enumerate(attributes):
                fields = entries.split(" - ")
                if fields[0] == 'COMMON-NAME':
                    common_name = fields[1].strip('')
                    enrxn_common_name.setdefault(unique_id, common_name)

    with open(reaction_path, 'r', encoding='iso-8859-1') as rp, open(metacyc_common_path, 'w',
                                                                     encoding='iso-8859-1') as mcp, open(
        common_name_not_found, 'w', encoding='iso-8859-1') as cnnf:
        mcp.write('# MetaCyc Version:\t' + version + '\n')
        for u_id, attrs in read_dat_itr(rp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            attributes = attrs.strip().split('\n')
            for index, entries in enumerate(attributes):
                fields = entries.split(" - ")
                if fields[0] == 'COMMON-NAME':
                    common_name = fields[1].strip('')
                    mcp.write(unique_id + '\t' + common_name + '\n')
                    break
                elif fields[0] == 'ENZYMATIC-REACTION':
                    enzymatic_reaction = fields[1].strip('')
                    try:
                        common_name = enrxn_common_name[enzymatic_reaction]
                        mcp.write(unique_id + '\t' + common_name + '\n')
                        break
                    except KeyError:
                        cnnf.write(u_id + '\n')
                        continue


def get_pathways_from_dat(file_path, version, output_folder):
    # input_folder = os.path.dirname(file_path)
    metacyc_pathways_name = "all_pwy.meta"
    metacyc_pathways_path = os.path.join(output_folder, metacyc_pathways_name)
    with open(file_path, 'r', encoding='iso-8859-1') as fp, open(metacyc_pathways_path, 'w',
                                                                 encoding='iso-8859-1') as mpp:
        mpp.write('# MetaCyc Version:\t' + version + '\n')
        for u_id, attrs in read_dat_itr(fp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            attributes = attrs.strip().split('\n')
            mpp.write(unique_id + '\n')


def get_ec_superseded_from_enzyme(file_path, version, output_folder):
    version = str(version)
    # input_folder = os.path.dirname(file_path)
    ec_super_name = 'EC-superseded'
    ec_super_path = os.path.join(output_folder, ec_super_name)
    ec_name_parsed_name = 'ec_name.map.parsed'
    ec_name_parsed_path = os.path.join(output_folder, ec_name_parsed_name)
    with open(file_path, 'r') as ip, open(ec_super_path, 'w') as esp, open(ec_name_parsed_path, 'w') as enp:
        esp.write('# ENZYME nomenclature database Version:\t' + version + '\n')
        enp.write('# ENZYME nomenclature database Version:\t' + version + '\n')
        for u_id, attrs in read_dat_itr(ip, "ID", " "):
            unique_id = 'EC-' + u_id.replace('ID', '').strip()
            attributes = attrs.strip().split('\n')
            common_name = []
            for entries in attributes:
                fields = entries.split()
                if fields[0] == 'DE' and 'Transferred entry:' in entries:
                    transferred_entry_line = entries.replace(' and ', ', ').split('Transferred entry:')
                    ec_nums = ['EC-' + ec_num.rstrip('.') for ec_num in transferred_entry_line[1].strip().split(', ')]
                    for ec_n in set(ec_nums):
                        esp.write(ec_n + '\t:TRANSFERRED-FROM\t' + unique_id + '\n')
                elif fields[0] == 'DE':
                    common_name.append(entries.replace('DE', '').strip())
            # Transfered EC Common-Name?
            if len(common_name) > 0:
                enp.write(u_id.replace('ID', '').strip() + '\t' + ''.join(common_name).rstrip('.') + '\n')


def find_default_version(pgdb_path):
    version_file = os.path.join(pgdb_path, 'default-version')
    if os.path.isfile(version_file):
        with open(version_file, 'r') as vf:
            for line in vf:
                version = line.strip()
                return version
    # print('default-version not found in', pgdb_path, ', exiting...')
    return None


def retrieve_exp_from_dat(file_path, exp_codes, taxon_list):
    exp_dict = {}
    with open(file_path, 'r', encoding='iso-8859-1') as fp:
        for u_id, attrs in read_dat_itr(fp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            # print(unique_id)
            attributes = attrs.strip().split('\n')
            found_exp = False
            found_taxon = False
            for index, entries in enumerate(attributes):
                fields = entries.split(" - ")
                if fields[0] == 'CITATIONS':
                    citations = fields[1].strip('')
                    exp_citations = [c for c in exp_codes if c in citations]
                    if len(exp_citations) > 0:
                        found_exp = True
                if fields[0] == 'SPECIES':
                    species_id = fields[1].strip('')
                    found_species_id = [s for s in taxon_list if s == species_id]
                    if len(found_species_id) > 0:
                        found_taxon = True
            if found_exp is True:
                if len(taxon_list) > 0 and found_taxon is True:
                    exp_dict.setdefault(unique_id, attrs)
                elif len(taxon_list) > 0 and found_taxon is False:
                    continue
                elif len(taxon_list) == 0:
                    exp_dict.setdefault(unique_id, attrs)

    return exp_dict


def refine_b_of_dat(reference_dat, input_dat, exp_codes, taxon_list):
    ref_exp_dict = retrieve_exp_from_dat(reference_dat, exp_codes, taxon_list)
    export_list = set()
    delete_list = set()
    with open(input_dat, 'r', encoding='iso-8859-1') as fp:
        input_uid = set()
        for u_id, attrs in read_dat_itr(fp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            input_uid.add(unique_id)
            # print(unique_id)
            attributes = attrs.strip().split('\n')
            found_exp = False
            for index, entries in enumerate(attributes):
                fields = entries.split(" - ")
                if fields[0] == 'CITATIONS':
                    citations = fields[1].strip('')
                    exp_citations = [c for c in exp_codes if c in citations]
                    if len(exp_citations) > 0:
                        found_exp = True
            if unique_id in ref_exp_dict.keys() and found_exp is False:
                exp_attrs = ref_exp_dict[unique_id]
                delete_list.add(unique_id)
                export_list.add(unique_id)

    export_list.update(set(ref_exp_dict.keys()) - input_uid)
    return sorted(export_list), sorted(delete_list)


def write_lisp_txt(pgdb_folder_name, id_list, lisp_task, input_type, output_dir):
    if lisp_task == 'export':
        with open(os.path.join(output_dir, pgdb_folder_name + '.' + input_type + '.export.txt'), 'w') as ep:
            ep.write('(in-package :ecocyc)\n')
            ep.write('(setq pwys \'(\n')
            for u_id in id_list:
                ep.write(u_id + '\n')
            ep.write('))\n')
    elif lisp_task == 'delete':
        with open(os.path.join(output_dir, pgdb_folder_name + '.' + input_type + '.delete.txt'), 'w') as ep:
            ep.write('(loop for p in \'(\n')
            for u_id in id_list:
                ep.write(u_id + '\n')
            ep.write(') do (delete-frame-and-dependents p))\n')


def read_protein_dat(file_path):
    enzrxn2species = {}
    with open(file_path, 'r', encoding='iso-8859-1') as fp:
        for u_id, attrs in read_dat_itr(fp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            attributes = attrs.strip().split('\n')
            enzrxn_set = set()
            species_set = set()

            for index, entries in enumerate(attributes):
                fields = entries.split(" - ")
                if fields[0] == 'CATALYZES':
                    enzrxn_id = fields[1].strip('')
                    enzrxn_set.add(enzrxn_id)
                if fields[0] == 'SPECIES':
                    species_id = fields[1].strip('')
                    species_set.add(species_id)
            if len(enzrxn_set) > 0 and len(species_set) > 0:
                for enzrxn_id in sorted(enzrxn_set):
                    try:
                        enzrxn2species[enzrxn_id] += list(species_set)
                    except KeyError:
                        enzrxn2species.setdefault(enzrxn_id, list(species_set))
    return enzrxn2species


def refine_enzrxn_list(enzrxn_list, enzrxn2species, taxon_list):
    refined_list = set()
    for enzrxn in enzrxn_list:
        try:
            species = enzrxn2species[enzrxn]
            if not set(species).isdisjoint(taxon_list):
                # print(enzrxn)
                refined_list.add(enzrxn)
        except KeyError:
            continue
    return sorted(refined_list)


def read_npp(file_path):
    npp_list = set()
    with open(file_path, 'r') as fp:
        for line in fp:
            if len(line) > 0 and line != '':
                npp_list.add(line.strip())
    return npp_list


def semi_auto_refine_b(reference_pgdb_folder, input_pgdb_folder, npp_file, taxon_list, exp_codes, output_dir):
    pgdb_name = os.path.basename(input_pgdb_folder)
    ref_pgdb_name = os.path.basename(reference_pgdb_folder)
    reference_pgdb_version = find_default_version(reference_pgdb_folder)
    input_pgdb_version = find_default_version(input_pgdb_folder)
    npp_list = read_npp(npp_file)
    if reference_pgdb_version is not None and input_pgdb_version is not None:
        input_pathways_dat_path = os.path.join(input_pgdb_folder, input_pgdb_version, 'data/pathways.dat')
        reference_pathways_dat_path = os.path.join(reference_pgdb_folder, reference_pgdb_version, 'data/pathways.dat')
        pathways_export_list, pathways_delete_list = refine_b_of_dat(reference_pathways_dat_path,
                                                                     input_pathways_dat_path,
                                                                     exp_codes, taxon_list)
        pathways_export_list = [p for p in pathways_export_list if p not in npp_list]
        if len(pathways_export_list) > 0:
            write_lisp_txt(pgdb_name, pathways_export_list, 'export', 'pathways', output_dir)
        if len(pathways_delete_list) > 0:
            write_lisp_txt(pgdb_name, pathways_delete_list, 'delete', 'pathways', output_dir)

        print('Total num. of pathways.dat EXP. UNIQUE-ID in ' + input_pgdb_folder + ' input to delete:',
              len(pathways_delete_list))
        print('Total num. of pathways.dat EXP. UNIQUE-ID in ' + reference_pgdb_folder + ' reference to export:',
              len(pathways_export_list))

        input_enzrxns_dat_path = os.path.join(input_pgdb_folder, input_pgdb_version, 'data/enzrxns.dat')
        reference_enzrxns_dat_path = os.path.join(reference_pgdb_folder, reference_pgdb_version, 'data/enzrxns.dat')
        enzrxns_export_list, enzrxn_delete_list = refine_b_of_dat(reference_enzrxns_dat_path, input_enzrxns_dat_path,
                                                                  exp_codes, [])
        if ref_pgdb_name != pgdb_name:
            reference_protein_dat_path = os.path.join(reference_pgdb_folder, reference_pgdb_version,
                                                      'data/proteins.dat')
            enzrxn2species = read_protein_dat(reference_protein_dat_path)
            enzrxns_export_list = refine_enzrxn_list(enzrxns_export_list, enzrxn2species, taxon_list)
            enzrxn_delete_list = refine_enzrxn_list(enzrxn_delete_list, enzrxn2species, taxon_list)

        if len(enzrxns_export_list) > 0:
            write_lisp_txt(pgdb_name, enzrxns_export_list, 'export', 'enzrxns', output_dir)
        if len(enzrxn_delete_list) > 0:
            write_lisp_txt(pgdb_name, enzrxn_delete_list, 'delete', 'enzrxns', output_dir)

        print('Total num. of enzrxns.dat EXP. UNIQUE-ID in ' + input_pgdb_folder + ' input to delete:',
              len(enzrxn_delete_list))
        print('Total num. of enzrxns.dat EXP. UNIQUE-ID in ' + reference_pgdb_folder + ' reference to export:',
              len(enzrxns_export_list))


def find_all_paths(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    paths = []
    if graph[start] not in path:
        new_paths = find_all_paths(graph, graph[start], end, path)
        for new_path in new_paths:
            paths.append(new_path)
    return paths


def get_children_taxons(classes_dat_path, taxon_id):
    graph = {}
    with open(classes_dat_path, 'r', encoding='iso-8859-1') as fp:
        for u_id, attrs in read_dat_itr(fp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            attributes = attrs.strip().split('\n')
            types = ""
            for index, entries in enumerate(attributes):
                fields = entries.split(" - ")
                if fields[0] == 'TYPES':
                    if fields[1].startswith('TAX-') or fields[1].startswith('ORG-'):
                        types = fields[1].strip('')
            if unique_id.startswith('ORG-') or unique_id.startswith('TAX-'):
                if types != "":
                    graph.setdefault(unique_id, types)
    if len(graph) > 0:
        childrens = []
        for cur_taxon in graph.keys():
            if cur_taxon not in childrens:
                paths_to_parents = find_all_paths(graph, cur_taxon, taxon_id)
                for paths in paths_to_parents:
                    childrens += paths
        return list(set(childrens))
    else:
        return []


def count_uniq_id_from_dat(dat_path):
    unique_id_set = set()
    with open(dat_path, 'r', encoding='iso-8859-1') as fp:
        for u_id, attrs in read_dat_itr(fp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            unique_id_set.add(unique_id)
    return len(unique_id_set)


def retrieve_polypeptides(file_path):
    polypeptides_list = set()
    with open(file_path, 'r', encoding='iso-8859-1') as fp:
        for u_id, attrs in read_dat_itr(fp, "UNIQUE-ID", " - "):
            unique_id = u_id.replace('UNIQUE-ID -', '').strip()
            attributes = attrs.strip().split('\n')
            for index, entries in enumerate(attributes):
                fields = entries.split(" - ")
                if fields[0].strip() == 'TYPES' and fields[1].strip() == 'Polypeptides':
                    polypeptides_list.add(unique_id)

    return sorted(polypeptides_list)


if __name__ == '__main__':
    parent_parser = ArgumentParser(description='Process PGDBs Flat Files', add_help=False)
    parent_parser.add_argument('-i', '--input', dest='input_path', required=True)
    parent_parser.add_argument('-o', '--output', dest='output_path', required=True)

    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest='task')
    subparsers.required = True
    reactions_parser = subparsers.add_parser('reactions', parents=[parent_parser])
    enzyme_parser = subparsers.add_parser('enzyme', parents=[parent_parser])
    metacyc_common_parser = subparsers.add_parser('common', parents=[parent_parser])
    metacyc_common_parser.add_argument('-e', '--enzrxn', dest='enzrxn_path', required=True)
    pathways_parser = subparsers.add_parser('pathways', parents=[parent_parser])
    refineb_parser = subparsers.add_parser('refineb', parents=[parent_parser])
    refineb_parser.add_argument('-r', '--reference', dest='ref_path', required=True)
    refineb_parser.add_argument('-t', '--tocreate', dest='to_create_path')
    refineb_parser.add_argument('-n', '--npp', dest='npp_path')
    stats_parser = subparsers.add_parser('stats', parents=[parent_parser])
    poly_parser = subparsers.add_parser('poly', parents=[parent_parser])
    # poly_parser.add_argument('-p', '--poly', dest='poly_path', required=True)

    args = parser.parse_args()

    input_path = args.input_path
    output_dir = args.output_path
    task_type = args.task

    if task_type == 'reactions':
        cur_version = get_versions_from_input(input_path)
        get_ec_mapping_from_reactions(input_path, cur_version, output_dir)
    elif task_type == 'enzyme':
        cur_version = get_versions_from_input(input_path)
        get_ec_superseded_from_enzyme(input_path, cur_version, output_dir)
    elif task_type == 'common':
        cur_version = get_versions_from_input(input_path)
        enzrxn_path = args.enzrxn_path
        get_metacyc_common_name(input_path, enzrxn_path, cur_version, output_dir)
    elif task_type == 'pathways':
        cur_version = get_versions_from_input(input_path)
        get_pathways_from_dat(input_path, cur_version, output_dir)
    elif task_type == 'refineb':
        ref_path = args.ref_path
        # output_dir = args.output_path
        exp_codes = ['EV-EXP', 'EV-AS', 'EV-IC']
        # print(input_path)
        # exp_dict = retrieve_exp_from_dat(input_path, exp_codes, [])
        # refine_b_of_dat(ref_path, input_path, exp_codes, ['TAX-4081'])
        # print(exp_dict.keys())
        if args.to_create_path is not None:
            # Metacyc to PGDB
            with open(args.to_create_path, 'r') as tcp:
                for line in tcp:
                    info = line.split('\t')
                    try:
                        pgdb_name = info[1].lower().strip()
                        if pgdb_name == 'aracyc':
                            taxon = 'ORG-' + info[5].strip()
                        else:
                            taxon = 'TAX-' + info[5].strip()
                        input_pgdb_folder = os.path.join(input_path, pgdb_name)
                        semi_auto_refine_b(ref_path, input_pgdb_folder, args.npp_file, [taxon], exp_codes, output_dir)
                    except IndexError:
                        continue
        else:
            ref_pgdb_name = os.path.basename(ref_path)
            input_pgdb_name = os.path.basename(input_path)
            if ref_pgdb_name == input_pgdb_name:
                # Old to New
                semi_auto_refine_b(ref_path, input_path, args.npp_file, [], exp_codes, output_dir)
            else:
                # Metacyc to Plantcyc
                reference_pgdb_version = find_default_version(ref_path)
                if reference_pgdb_version is not None:
                    reference_classes_dat_path = os.path.join(ref_path, reference_pgdb_version, 'data/classes.dat')
                    plant_taxons = get_children_taxons(reference_classes_dat_path, "TAX-33090")
                    print('Found Children Taxon Num.:\t', len(plant_taxons))
                    semi_auto_refine_b(ref_path, input_path, args.npp_file, plant_taxons, exp_codes, output_dir)

    elif task_type == 'stats':
        # generate stats of pgdbs
        input_pgdb_user_path = input_path
        out_put_path = os.path.join(output_dir, 'pgdbs_stats.txt')
        with open(out_put_path, 'w') as op:
            op.write('PGDB\tPathways\tEnzymes\tReactions\tCompounds\n')
            for cyc in sorted(os.listdir(input_path)):
                cyc_path = os.path.join(input_path, cyc)
                if not os.path.isdir(cyc_path):
                    continue
                default_version = find_default_version(cyc_path)
                if default_version is not None:
                    print('Processing PGDB:', cyc)
                    pathways_dat_path = os.path.join(cyc_path, default_version, 'data/pathways.dat')
                    proteins_dat_path = os.path.join(cyc_path, default_version, 'data/proteins.dat')
                    reactions_dat_path = os.path.join(cyc_path, default_version, 'data/reactions.dat')
                    compounds_dat_path = os.path.join(cyc_path, default_version, 'data/compounds.dat')
                    if os.path.isfile(pathways_dat_path):
                        num_of_pathways = str(count_uniq_id_from_dat(pathways_dat_path))
                    else:
                        num_of_pathways = 'NA'
                    if os.path.isfile(proteins_dat_path):
                        num_of_proteins = str(count_uniq_id_from_dat(proteins_dat_path))
                    else:
                        num_of_proteins = 'NA'
                    if os.path.isfile(reactions_dat_path):
                        num_of_reactions = str(count_uniq_id_from_dat(reactions_dat_path))
                    else:
                        num_of_reactions = 'NA'
                    if os.path.isfile(compounds_dat_path):
                        num_of_compounds = str(count_uniq_id_from_dat(compounds_dat_path))
                    else:
                        num_of_compounds = 'NA'
                else:
                    num_of_pathways = 'NA'
                    num_of_proteins = 'NA'
                    num_of_reactions = 'NA'
                    num_of_compounds = 'NA'

                cyc_name = re.sub('[Cc]yc$', 'Cyc', cyc.title())
                op.write('\t'.join(
                    [cyc_name + ' ' + default_version, num_of_pathways, num_of_proteins, num_of_reactions,
                     num_of_compounds]) + '\n')
    elif task_type == 'poly':
        input_pgdb_protein_dat_path = input_path
        out_put_path = os.path.join(output_dir, 'polypeptides_list.txt')
        polypeptides_list = retrieve_polypeptides(input_pgdb_protein_dat_path)
        with open(out_put_path, 'w') as op:
            for poly in polypeptides_list:
                op.write(poly + '\n')
