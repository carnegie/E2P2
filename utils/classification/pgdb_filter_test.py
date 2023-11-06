import os
import re

from Bio import SeqIO
from tqdm import tqdm


def read_biocyc_flat_file(fp):
    unique_id, attrs = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith("UNIQUE-ID"):
            if unique_id:
                yield unique_id, attrs
            unique_id, attrs = re.sub(r'^UNIQUE-ID[\s|\t]+-[\s|\t]+', '', line), []
        else:
            info = [i.strip() for i in re.split(r'[\s|\t]+-[\s|\t]+', line)]
            try:
                if 'IN-PATHWAY' in info[0] or 'REACTION-LIST' in info[0]:
                    attrs.append(info[1])
            except IndexError:
                continue
    if unique_id:
        yield unique_id, attrs


def process_reactions_dat(reactions_dat_path):
    reaction_to_pathways = {}
    with open(reactions_dat_path, 'r', encoding='iso-8859-1') as rdp:
        for unique_id, pathways in read_biocyc_flat_file(rdp):
            try:
                reaction_to_pathways[unique_id].update(pathways)
            except KeyError:
                reaction_to_pathways.setdefault(unique_id, set(pathways))
    return reaction_to_pathways


def process_pathways_dat(pathways_dat_path):
    reaction_to_pathways = {}
    with open(pathways_dat_path, 'r', encoding='iso-8859-1') as rdp:
        for unique_id, reactions in read_biocyc_flat_file(rdp):
            for reaction in sorted(reactions):
                try:
                    reaction_to_pathways[reaction].add(unique_id)
                except KeyError:
                    reaction_to_pathways.setdefault(reaction, {unique_id})
    return reaction_to_pathways


def read_pf_file(fp):
    seq_id, metacyc_rxns = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith("ID"):
            if seq_id:
                yield seq_id, metacyc_rxns
            seq_id, metacyc_rxns = re.sub(r'^ID\t', '', line), []
        else:
            info = [i.strip() for i in re.split(r'\t', line)]
            try:
                if 'METACYC' in info[0]:
                    metacyc_rxns.append(info[1])
            except IndexError:
                continue
    if seq_id:
        yield seq_id, metacyc_rxns


def process_pf_file(pf_file_path):
    rxn_ids = {}
    with open(pf_file_path, 'r') as pfp:
        for uniq_id, attrs in read_pf_file(pfp):
            for attr in attrs:
                try:
                    rxn_ids[attr].add(uniq_id)
                except KeyError:
                    rxn_ids.setdefault(attr, {uniq_id})
    return rxn_ids


def process_rpsd(rpsd_path, ec_superseded_path=None):
    ec_superseded = {}
    if ec_superseded_path:
        with open(ec_superseded_path, 'r') as esp:
            for line in esp:
                if line.startswith('#'):
                    continue
                try:
                    superseded_ec = re.sub('EC-', '', re.search(r'^EC-(\d+\.){3}\w+', line).group(0))
                    original_ec = re.sub('EC-', '', re.search(r'EC-(\d+\.){3}\w+$', line).group(0))
                    try:
                        ec_superseded[original_ec].add(superseded_ec)
                    except KeyError:
                        ec_superseded.setdefault(original_ec, {superseded_ec})
                except AttributeError:
                    print('Line Error', AttributeError)
    seq_to_annot = {}
    annot_to_seq = {}
    with open(rpsd_path, 'r') as rp:
        for record in SeqIO.parse(rp, 'fasta'):
            header = record.id
            info = [i.strip() for i in header.split('|')]
            uniprot_id = info[0]
            annotations = [a for a in info[1:] if a != '']
            annotations_updated = set()
            for annot in annotations:
                try:
                    annotations_updated.update(ec_superseded[annot])
                except KeyError:
                    annotations_updated.add(annot)
            try:
                seq_to_annot[uniprot_id].update(set(annotations_updated))
            except KeyError:
                seq_to_annot.setdefault(uniprot_id, set(annotations_updated))
            for annot in sorted(annotations_updated):
                try:
                    annot_to_seq[annot].add(uniprot_id)
                except KeyError:
                    annot_to_seq.setdefault(annot, {uniprot_id})
    return seq_to_annot, annot_to_seq


def read_ef_mapping(ef_mapping_path, ec_superseded_path=None):
    ec_superseded = {}
    if ec_superseded_path:
        with open(ec_superseded_path, 'r') as esp:
            for line in esp:
                if line.startswith('#'):
                    continue
                try:
                    superseded_ec = re.sub('EC-', '', re.search(r'^EC-(\d+\.){3}\w+', line).group(0))
                    original_ec = re.sub('EC-', '', re.search(r'EC-(\d+\.){3}\w+$', line).group(0))
                    try:
                        ec_superseded[original_ec].add(superseded_ec)
                    except KeyError:
                        ec_superseded.setdefault(original_ec, {superseded_ec})
                except AttributeError:
                    print('Line Error', AttributeError)

    ef_to_id = {}
    id_to_ef = {}
    with open(ef_mapping_path, 'r') as emp:
        for line in emp:
            if line.startswith('#'):
                continue
            info = [i.strip() for i in line.split('\t')]
            try:
                efclass = info[0]
                fun_id = [f for f in info[1].split('|') if f != '']
                fun_id_updated = set()
                for f in fun_id:
                    try:
                        fun_id_updated.update(ec_superseded[f])
                    except KeyError:
                        fun_id_updated.add(f)
                ef_to_id.setdefault(efclass, set(fun_id_updated))
                for f in sorted(fun_id_updated):
                    try:
                        id_to_ef[f].add(efclass)
                    except KeyError:
                        id_to_ef.setdefault(f, {efclass})
            except IndexError:
                print('Line parsing error', line)
                continue
    return ef_to_id, id_to_ef


def filter_performances(performance_path, threshold=0.5, score_index=-1):
    kept_classes = set()
    with open(performance_path, 'r') as pp:
        for line in pp:
            info = [i.strip() for i in line.split('\t')]
            try:
                input_class = info[0]
                score = float(info[score_index])
                if score >= threshold:
                    kept_classes.add(input_class)
            except (IndexError, ValueError) as e:
                print('Line parsing error', line)
    return kept_classes


def partition_annotations(annot_to_seq):
    ec_1_to_seq = {}
    ec_2_to_seq = {}
    ec_3_to_seq = {}
    ec_to_seq = {}
    rxn_to_seq = {}
    for annot in annot_to_seq:
        seqs_of_annot = set(annot_to_seq[annot])
        if re.search(r'^(\d+\.){3}\d+$', annot):
            ec = annot.split('.')
            if len(ec) != 4:
                print('EC Error', annot)
                continue
            try:
                ec_1_to_seq[ec[0]].update(set(seqs_of_annot))
            except KeyError:
                ec_1_to_seq.setdefault(ec[0], set(seqs_of_annot))
            try:
                ec_2_to_seq['.'.join(ec[:2])].update(set(seqs_of_annot))
            except KeyError:
                ec_2_to_seq.setdefault('.'.join(ec[:2]), set(seqs_of_annot))
            try:
                ec_3_to_seq['.'.join(ec[:3])].update(set(seqs_of_annot))
            except KeyError:
                ec_3_to_seq.setdefault('.'.join(ec[:3]), set(seqs_of_annot))
            ec_to_seq.setdefault(annot, set(seqs_of_annot))
        else:
            rxn_to_seq.setdefault(annot, set(seqs_of_annot))
    return ec_1_to_seq, ec_2_to_seq, ec_3_to_seq, ec_to_seq, rxn_to_seq


def read_metacyc_rxn_to_ec_maps(official_ec_to_rxn_path, rxn_to_ec_path, ec_list=None, rxn_list=None):
    official_rxn_to_ec_dict = {}
    others_rxn_to_ec_dict = {}
    with open(official_ec_to_rxn_path, 'r') as oetrp:
        for line in oetrp:
            if line.startswith('#'):
                continue
            info = [i.strip() for i in line.split('\t')]
            try:
                ec_num = info[0]
                rxn_id = info[1]
                if ec_list and ec_num not in ec_list:
                    continue
                if rxn_list and rxn_id not in rxn_list:
                    continue
                if not re.search(r'^(\d+\.){3}\d+$', ec_num):
                    continue
                try:
                    official_rxn_to_ec_dict[rxn_id].add(ec_num)
                except KeyError:
                    official_rxn_to_ec_dict.setdefault(rxn_id, {ec_num})
            except IndexError:
                print('Line parsing error', info)
                continue
    with open(rxn_to_ec_path, 'r') as rtep:
        for line in rtep:
            if line.startswith('#'):
                continue
            info = [i.strip() for i in line.split('\t')]
            try:
                ec_num = info[1]
                rxn_id = info[0]
                if ec_list and ec_num not in ec_list:
                    continue
                if rxn_list and rxn_id not in rxn_list:
                    continue
                if not re.search(r'^(\d+\.){3}\d+$', ec_num):
                    continue
                if rxn_id in official_rxn_to_ec_dict:
                    if ec_num in official_rxn_to_ec_dict[rxn_id]:
                        continue
                try:
                    others_rxn_to_ec_dict[rxn_id].add(ec_num)
                except KeyError:
                    others_rxn_to_ec_dict.setdefault(rxn_id, {ec_num})
            except IndexError:
                print('Line parsing error', info)
                continue
    return official_rxn_to_ec_dict, others_rxn_to_ec_dict


def filter_jaccard_rpsd(jaccard_kept_classes, jaccard_ef_to_id, annot_to_seq, official_rxn_to_ec_dict,
                        others_rxn_to_ec_dict):
    jaccard_kept_ids = set()
    for ef in sorted(jaccard_kept_classes):
        try:
            mapped_ids = jaccard_ef_to_id[ef]
            jaccard_kept_ids.update(mapped_ids)
        except KeyError:
            print('EF not in map', ef)
    jaccard_kept_ec_seqs = set()
    jaccard_kept_ec_to_seq = {}
    jaccard_kept_rxn_seqs = set()
    jaccard_kept_rxn_to_seq = {}
    for annot in sorted(jaccard_kept_ids):
        try:
            seq_of_annot = annot_to_seq[annot]
            if re.search(r'^(\d+\.){3}\d+$', annot):
                jaccard_kept_ec_seqs.update(seq_of_annot)
                jaccard_kept_ec_to_seq.setdefault(annot, set(seq_of_annot))
            else:
                jaccard_kept_rxn_seqs.update(seq_of_annot)
                jaccard_kept_rxn_to_seq.setdefault(annot, set(seq_of_annot))
        except KeyError:
            print('No seq found', annot)
    print('Jaccard Class', len(jaccard_kept_classes), 'Jaccard IDs', len(jaccard_kept_ids),
          'Jaccard Seqs', len(jaccard_kept_ec_seqs | jaccard_kept_rxn_seqs))
    print('Jaccard RXN', len(jaccard_kept_rxn_to_seq), 'Jaccard RXN Seqs',
          len(jaccard_kept_rxn_seqs), 'Jaccard EC', len(jaccard_kept_ec_to_seq), 'Jaccard EC Seqs',
          len(jaccard_kept_ec_seqs))
    for ec in sorted(jaccard_kept_ec_to_seq.keys()):
        for rxn in official_rxn_to_ec_dict:
            if ec in official_rxn_to_ec_dict[rxn]:
                jaccard_kept_rxn_seqs.update(set(jaccard_kept_ec_to_seq[ec]))
                try:
                    jaccard_kept_rxn_to_seq[rxn].update(set(jaccard_kept_ec_to_seq[ec]))
                except KeyError:
                    jaccard_kept_rxn_to_seq.setdefault(rxn, set(jaccard_kept_ec_to_seq[ec]))
        for rxn in others_rxn_to_ec_dict:
            if ec in others_rxn_to_ec_dict[rxn]:
                jaccard_kept_rxn_seqs.update(set(jaccard_kept_ec_to_seq[ec]))
                try:
                    jaccard_kept_rxn_to_seq[rxn].update(set(jaccard_kept_ec_to_seq[ec]))
                except KeyError:
                    jaccard_kept_rxn_to_seq.setdefault(rxn, set(jaccard_kept_ec_to_seq[ec]))
    print('Jaccard RXN w/ EC mapped', len(jaccard_kept_rxn_to_seq),
          'Jaccard RXN w/ EC mapped Seqs', len(jaccard_kept_rxn_seqs))
    return jaccard_kept_rxn_to_seq


def filter_main(prefix):
    size = 10
    threshold = 0.5
    rpsd_path = '/Users/bxue/Documents/Carnegie/NLP/Test/fasta/rpsd-4.2-20190307.nlp.fasta'
    ec_superseded_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/release_2019-03-07/maps/' \
                         'pf-EC-superseded.mapping'
    official_ec_to_rxn_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/release_2019-03-07/maps/' \
                              'pf-official-EC-metacyc-RXN.mapping'
    rxn_to_ec_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/release_2019-03-07/maps/' \
                     'pf-metacyc-RXN-EC.mapping'

    jaccard_rpsd = '/Users/bxue/Documents/Carnegie/NLP/Test/fasta/rpsd-4.2-20190307.nlp.jaccard.fasta'
    jaccard_ef_mapping = '//Users/bxue/Documents/Carnegie/NLP/Test/fasta/efclasses.nlp.mapping'
    jaccard_performance = '/Users/bxue/Documents/Carnegie/NLP/Test/fasta/performance/' \
                          'rpsd-4.2-20190307.nlp.jaccard.label.performance_individual.with_size.txt'

    seq_to_annot, annot_to_seq = process_rpsd(rpsd_path, ec_superseded_path)
    print('RPSD Seq', len(seq_to_annot), 'RPSD Annot', len(annot_to_seq))

    official_rxn_to_ec_dict, others_rxn_to_ec_dict = \
        read_metacyc_rxn_to_ec_maps(official_ec_to_rxn_path, rxn_to_ec_path)

    ec_1_to_seq, ec_2_to_seq, ec_3_to_seq, ec_to_seq, rxn_to_seq = partition_annotations(annot_to_seq)
    print('RPSD EC', len(ec_to_seq), 'RPSD RXN', len(rxn_to_seq))
    rxn_with_ec_mapped_to_seq = {}

    all_rxn_seqs = set()
    over_size_rxn_seqs = set()
    for rxn in sorted(rxn_to_seq.keys()):
        rxn_with_ec_mapped_to_seq.setdefault(rxn, rxn_to_seq[rxn])
        all_rxn_seqs.update(set(rxn_to_seq[rxn]))
        if len(rxn_to_seq[rxn]) >= size:
            over_size_rxn_seqs.update(set(rxn_to_seq[rxn]))
        else:
            rxn_to_seq.pop(rxn, None)
    print('All RXN Seq', len(all_rxn_seqs), '; >=', size, 'RXN Seq', len(over_size_rxn_seqs), '; >=', size, 'RXN',
          len(rxn_to_seq))

    all_ec_seqs = set()
    over_size_ec_seqs = set()
    for ec in sorted(ec_to_seq.keys()):
        all_ec_seqs.update(set(ec_to_seq[ec]))
        for rxn in official_rxn_to_ec_dict:
            if ec in official_rxn_to_ec_dict[rxn]:
                try:
                    rxn_with_ec_mapped_to_seq[rxn].update(set(ec_to_seq[ec]))
                except KeyError:
                    rxn_with_ec_mapped_to_seq.setdefault(rxn, set(ec_to_seq[ec]))
        for rxn in others_rxn_to_ec_dict:
            if ec in others_rxn_to_ec_dict[rxn]:
                try:
                    rxn_with_ec_mapped_to_seq[rxn].update(set(ec_to_seq[ec]))
                except KeyError:
                    rxn_with_ec_mapped_to_seq.setdefault(rxn, set(ec_to_seq[ec]))
        if len(set(ec_to_seq[ec])) >= size:
            over_size_ec_seqs.update(set(ec_to_seq[ec]))
        else:
            ec_to_seq.pop(ec, None)

    print('All EC Seq', len(all_ec_seqs), '; >=', size, 'EC Seq', len(over_size_ec_seqs), '; >=', size, 'EC',
          len(ec_to_seq))
    print('All >=', size, 'Seq', len(over_size_ec_seqs | over_size_rxn_seqs))
    print('RPSD RXN w/ EC mapped', len(rxn_with_ec_mapped_to_seq))
    all_rxn_with_ec_mapped_seqs = set()
    over_size_all_rxn_with_ec_mapped_seqs = set()
    for rxn in sorted(rxn_with_ec_mapped_to_seq.keys()):
        all_rxn_with_ec_mapped_seqs.update(set(rxn_with_ec_mapped_to_seq[rxn]))
        if len(set(rxn_with_ec_mapped_to_seq[rxn])) >= size:
            over_size_all_rxn_with_ec_mapped_seqs.update(set(rxn_with_ec_mapped_to_seq[rxn]))
        else:
            rxn_with_ec_mapped_to_seq.pop(rxn, None)
    print('All RXN w/ EC mapped Seq', len(all_rxn_with_ec_mapped_seqs), '; >=', size, 'RXN w/ EC mapped Seq',
          len(over_size_all_rxn_with_ec_mapped_seqs), '; >=', size, 'RXN w/ EC mapped', len(rxn_with_ec_mapped_to_seq))

    jaccard_ef_to_id, jacard_id_to_ef = read_ef_mapping(jaccard_ef_mapping, ec_superseded_path)
    jaccard_seq_to_annot, jaccard_annot_to_seq = process_rpsd(jaccard_rpsd)
    over_size_jaccard_kept_classes = set()
    over_size_jaccard_kept_seqs = set()
    for annot in jaccard_annot_to_seq:
        if len(jaccard_annot_to_seq[annot]) >= size:
            over_size_jaccard_kept_classes.add(annot)
            over_size_jaccard_kept_seqs.update(jaccard_annot_to_seq[annot])

    print('Filter Jaccard by Size >=', size)
    over_size_jaccard_kept_rxn_to_seq = filter_jaccard_rpsd(over_size_jaccard_kept_classes, jaccard_ef_to_id,
                                                                       annot_to_seq, official_rxn_to_ec_dict,
                                                                       others_rxn_to_ec_dict)
    if prefix == 'Fscore':
        performance_jaccard_kept_classes = filter_performances(jaccard_performance, threshold=threshold, score_index=4)
        print('Fscore', len(performance_jaccard_kept_classes))
    elif prefix == 'Precision':
        performance_jaccard_kept_classes = filter_performances(jaccard_performance, threshold=threshold, score_index=2)
        print('Precision', len(performance_jaccard_kept_classes))
    else:
        print('Prefix error', prefix)
        raise SystemError
    print('Filter Jaccard by Performance >=', threshold)
    performance_jaccard_kept_rxn_to_seq = filter_jaccard_rpsd(performance_jaccard_kept_classes, jaccard_ef_to_id,
                                                              annot_to_seq, official_rxn_to_ec_dict,
                                                              others_rxn_to_ec_dict)

    return rxn_with_ec_mapped_to_seq, over_size_jaccard_kept_rxn_to_seq, performance_jaccard_kept_rxn_to_seq


def process_pgdb_folder(pgdb_folder_path, plantcyc=False):
    reactions_dat_reaction_to_pathway = set()
    pathways_dat_reaction_to_pathway = set()
    pf_reactions = set()
    for f in os.listdir(pgdb_folder_path):
        file_path = os.path.join(pgdb_folder_path, f)
        if os.path.isdir(file_path):
            reactions_dat_path = os.path.join(file_path, 'data', 'reactions.dat')
            pathways_dat_path = os.path.join(file_path, 'data', 'pathways.dat')

            pf_file_path = os.path.join(file_path, 'input')
            if os.path.isdir(pf_file_path):
                for pf in os.listdir(pf_file_path):
                    if pf.endswith('.pf'):
                        pf_file_path = os.path.join(pf_file_path, pf)
                        break
            if os.path.isfile(reactions_dat_path) and os.path.isfile(pathways_dat_path):
                if not plantcyc and os.path.isfile(pf_file_path):
                    reactions_dat_reaction_to_pathway = process_reactions_dat(reactions_dat_path)
                    pathways_dat_reaction_to_pathway = process_pathways_dat(pathways_dat_path)
                    pf_reactions = process_pf_file(pf_file_path)
                elif plantcyc:
                    reactions_dat_reaction_to_pathway = process_reactions_dat(reactions_dat_path)
                    pathways_dat_reaction_to_pathway = process_pathways_dat(pathways_dat_path)
                else:
                    print('File error', pf_file_path, file_path)
                    raise SystemError
            else:
                print('Folder error', file_path)
                raise SystemError
    return reactions_dat_reaction_to_pathway, pathways_dat_reaction_to_pathway, pf_reactions


def key_to_set_dict_values(key_set_dict):
    values = set()
    for key in key_set_dict:
        values.update(key_set_dict[key])
    return values


def combine_list_of_key_to_set_dicts(list_of_dicts, output_dict=None):
    if output_dict is None:
        output_dict = {}
    for d in list_of_dicts:
        for key in d:
            try:
                output_dict[key].update(set(d[key]))
            except KeyError:
                output_dict.setdefault(key, set(d[key]))
    return output_dict


def pgdb_helper_1(prefix, rxn_list, reactions_dat_reaction_to_pathway, pathways_dat_reaction_to_pathway, pf_reactions,
                  all_reactions_dat_reaction_to_pathway, all_pathways_dat_reaction_to_pathway, all_pf_reactions):
    combine_list_of_key_to_set_dicts([reactions_dat_reaction_to_pathway],
                                     all_reactions_dat_reaction_to_pathway)
    combine_list_of_key_to_set_dicts([pathways_dat_reaction_to_pathway],
                                     all_pathways_dat_reaction_to_pathway)
    all_pf_reactions.update(pf_reactions)

    reactions_dat_reactions = set(reactions_dat_reaction_to_pathway.keys())
    pathways_dat_reactions = set(pathways_dat_reaction_to_pathway.keys())
    reactions_dat_pathways = key_to_set_dict_values(reactions_dat_reaction_to_pathway)
    pathways_dat_pathways = key_to_set_dict_values(pathways_dat_reaction_to_pathway)

    reactions_dat_reaction_to_pathway_subset = {k: reactions_dat_reaction_to_pathway[k] for k in
                                                sorted(set(reactions_dat_reaction_to_pathway.keys()) & set(rxn_list))}
    pathways_dat_reaction_to_pathway_subset = {k: pathways_dat_reaction_to_pathway[k] for k in
                                               sorted(set(pathways_dat_reaction_to_pathway.keys()) & set(rxn_list))}
    reactions_dat_reactions_subset = set(reactions_dat_reaction_to_pathway_subset.keys())
    pathways_dat_reactions_subset = set(pathways_dat_reaction_to_pathway_subset.keys())
    reactions_dat_pathways_subset = key_to_set_dict_values(reactions_dat_reaction_to_pathway_subset)
    pathways_dat_pathways_subset = key_to_set_dict_values(pathways_dat_reaction_to_pathway_subset)

    pf_reactions_subset = set(k for k in sorted(pf_reactions) if k in rxn_list)

    print(prefix, len(reactions_dat_reactions), len(reactions_dat_reactions_subset),
          len(pathways_dat_reactions), len(pathways_dat_reactions_subset),
          len(reactions_dat_pathways), len(reactions_dat_pathways_subset),
          len(pathways_dat_pathways), len(pathways_dat_pathways_subset), len(pf_reactions), len(pf_reactions_subset))
    # return \
    #     reactions_dat_reactions, pathways_dat_reactions, reactions_dat_pathways, pathways_dat_pathways, \
    #     reactions_dat_reactions_subset, pathways_dat_reactions_subset, reactions_dat_pathways_subset, \
    #     pathways_dat_pathways_subset, pf_reactions_subset
    return reactions_dat_reactions, reactions_dat_reactions_subset, pathways_dat_pathways, pathways_dat_pathways_subset, \
           reactions_dat_reaction_to_pathway_subset, pathways_dat_reaction_to_pathway_subset


def pgdb_helper_2(prefix, rxn_list, reactions_dat_reaction_to_pathway, pathways_dat_reaction_to_pathway, pf_reactions):
    reactions_dat_reactions = set(reactions_dat_reaction_to_pathway.keys())
    pathways_dat_reactions = set(pathways_dat_reaction_to_pathway.keys())
    reactions_dat_pathways = key_to_set_dict_values(reactions_dat_reaction_to_pathway)
    pathways_dat_pathways = key_to_set_dict_values(pathways_dat_reaction_to_pathway)

    reactions_dat_reaction_to_pathway_subset = {k: reactions_dat_reaction_to_pathway[k] for k
                                                             in
                                                             sorted(set(
                                                                 reactions_dat_reaction_to_pathway.keys()) & set(
                                                                 rxn_list))}
    pathways_dat_reaction_to_pathway_subset = {k: pathways_dat_reaction_to_pathway[k] for k in
                                                            sorted(set(
                                                                pathways_dat_reaction_to_pathway.keys()) & set(
                                                                rxn_list))}
    reactions_dat_reactions_subset = set(reactions_dat_reaction_to_pathway_subset.keys())
    pathways_dat_reactions_subset = set(pathways_dat_reaction_to_pathway_subset.keys())
    reactions_dat_pathways_subset = key_to_set_dict_values(
        reactions_dat_reaction_to_pathway_subset)
    pathways_dat_pathways_subset = key_to_set_dict_values(
        pathways_dat_reaction_to_pathway_subset)
    pf_reactions_subset = set(k for k in sorted(pf_reactions) if k in rxn_list)
    print(prefix, len(reactions_dat_reactions),
          len(reactions_dat_reactions_subset),
          len(pathways_dat_reactions), len(pathways_dat_reactions_subset),
          len(reactions_dat_pathways), len(reactions_dat_pathways_subset),
          len(pathways_dat_pathways), len(pathways_dat_pathways_subset),
          len(pf_reactions), len(pf_reactions_subset))


def pgdb_main(prefix, rxn_list):
    output_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/filtering_wo_savi'
    plantcyc_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/plantcyc_02_16_21'
    pre_savi_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_preSAVI_02_16_21'
    release_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_release_02_16_21'

    plantcyc_reactions_dat_reaction_to_pathway, plantcyc_pathways_dat_reaction_to_pathway, _ = \
        process_pgdb_folder(plantcyc_folder, plantcyc=True)
    plantcyc_reactions_dat_reactions = set(plantcyc_reactions_dat_reaction_to_pathway.keys())
    plantcyc_pathways_dat_reactions = set(plantcyc_pathways_dat_reaction_to_pathway.keys())
    plantcyc_reactions_dat_pathways = key_to_set_dict_values(plantcyc_reactions_dat_reaction_to_pathway)
    plantcyc_pathways_dat_pathways = key_to_set_dict_values(plantcyc_pathways_dat_reaction_to_pathway)

    plantcyc_reactions_dat_reaction_to_pathway_subset = {k: plantcyc_reactions_dat_reaction_to_pathway[k] for k in
                                                         sorted(set(
                                                             plantcyc_reactions_dat_reaction_to_pathway.keys()) & set(
                                                             rxn_list))}
    plantcyc_pathways_dat_reaction_to_pathway_subset = {k: plantcyc_pathways_dat_reaction_to_pathway[k] for k in
                                                        sorted(
                                                            set(plantcyc_pathways_dat_reaction_to_pathway.keys()) & set(
                                                                rxn_list))}
    plantcyc_reactions_dat_reactions_subset = set(plantcyc_reactions_dat_reaction_to_pathway_subset.keys())
    plantcyc_pathways_dat_reactions_subset = set(plantcyc_pathways_dat_reaction_to_pathway_subset.keys())
    plantcyc_reactions_dat_pathways_subset = key_to_set_dict_values(plantcyc_reactions_dat_reaction_to_pathway_subset)
    plantcyc_pathways_dat_pathways_subset = key_to_set_dict_values(plantcyc_pathways_dat_reaction_to_pathway_subset)

    print(prefix, 'plantcyc', len(plantcyc_reactions_dat_reactions), len(plantcyc_reactions_dat_reactions_subset),
          len(plantcyc_pathways_dat_reactions), len(plantcyc_pathways_dat_reactions_subset),
          len(plantcyc_reactions_dat_pathways), len(plantcyc_reactions_dat_pathways_subset),
          len(plantcyc_pathways_dat_pathways), len(plantcyc_pathways_dat_pathways_subset))

    all_pre_savi_reactions_dat_reaction_to_pathway, all_pre_savi_pathways_dat_reaction_to_pathway, all_pre_savi_pf_reactions = \
        {}, {}, set()
    all_release_reactions_dat_reaction_to_pathway, all_release_pathways_dat_reaction_to_pathway, all_release_pf_reactions = \
        {}, {}, set()
    for f in sorted(set(os.listdir(pre_savi_folder)) & set(os.listdir(release_folder))):
        if f.startswith('.'):
            continue
        # if not f.startswith('aracyc') and not f.startswith('oryzacyc'):
        #     continue
        pre_savi_pgdb = os.path.join(pre_savi_folder, f)
        pre_savi_reactions_dat_reaction_to_pathway, pre_savi_pathways_dat_reaction_to_pathway, pre_savi_pf_reactions = \
            process_pgdb_folder(pre_savi_pgdb)

        pre_savi_reactions_dat_reactions, pre_savi_reactions_dat_reactions_subset, pre_savi_pathways_dat_pathways, \
        pre_savi_pathways_dat_pathways_subset, pre_savi_reactions_dat_reaction_to_pathway_subset, pre_savi_pathways_dat_reaction_to_pathway_subset = \
            pgdb_helper_1(prefix + ' PreSAVI ' + f + '_PreSAVI', rxn_list, pre_savi_reactions_dat_reaction_to_pathway,
                          pre_savi_pathways_dat_reaction_to_pathway, pre_savi_pf_reactions,
                          all_pre_savi_reactions_dat_reaction_to_pathway, all_pre_savi_pathways_dat_reaction_to_pathway,
                          all_pre_savi_pf_reactions)

        release_pgdb = os.path.join(release_folder, f)
        release_reactions_dat_reaction_to_pathway, release_pathways_dat_reaction_to_pathway, release_pf_reactions = \
            process_pgdb_folder(release_pgdb)
        release_reactions_dat_reactions, release_reactions_dat_reactions_subset, release_pathways_dat_pathways, \
        release_pathways_dat_pathways_subset, release_reactions_dat_reaction_to_pathway_subset, release_pathways_dat_reaction_to_pathway_subset = \
            pgdb_helper_1(prefix + ' Release ' + f + '_Release', rxn_list, release_reactions_dat_reaction_to_pathway,
                          release_pathways_dat_reaction_to_pathway, release_pf_reactions,
                          all_release_reactions_dat_reaction_to_pathway, all_release_pathways_dat_reaction_to_pathway,
                          all_release_pf_reactions)

        yield f, pre_savi_reactions_dat_reaction_to_pathway, pre_savi_pathways_dat_reaction_to_pathway, pre_savi_reactions_dat_reaction_to_pathway_subset, pre_savi_pathways_dat_reaction_to_pathway_subset, \
              release_reactions_dat_reaction_to_pathway, release_reactions_dat_reaction_to_pathway_subset, release_pathways_dat_reaction_to_pathway_subset, release_pathways_dat_reaction_to_pathway
        # with open(os.path.join(output_folder, re.sub(r'[><=]+', '_', prefix) + '.PreSAVI.' + f + '.reaction.txt'), 'w') as op1:
        #     op1.write('RXN\tKept\n')
        #     for rxn in sorted(pre_savi_reactions_dat_reactions | pre_savi_reactions_dat_reactions_subset):
        #         if rxn in pre_savi_reactions_dat_reactions_subset:
        #             op1.write(rxn + '\tY\n' )
        #         else:
        #             op1.write(rxn + '\tN\n')
        # with open(os.path.join(output_folder, re.sub(r'[><=]+', '_', prefix) + '.PreSAVI.' + f + '.pathway.txt'), 'w') as op2:
        #     op2.write('Pathway\tKept\n')
        #     for pwy in sorted(pre_savi_pathways_dat_pathways | pre_savi_pathways_dat_pathways_subset):
        #         if pwy in pre_savi_pathways_dat_pathways_subset:
        #             op2.write(pwy + '\tY\n' )
        #         else:
        #             op2.write(pwy + '\tN\n')

        # with open(os.path.join(output_folder, re.sub(r'[><=]+', '_', prefix) + '.Release.' + f + '.reaction.txt'), 'w') as op1:
        #     op1.write('RXN\tKept\n')
        #     for rxn in sorted(release_reactions_dat_reactions | release_reactions_dat_reactions_subset):
        #         if rxn in release_reactions_dat_reactions_subset:
        #             op1.write(rxn + '\tY\n' )
        #         else:
        #             op1.write(rxn + '\tN\n')
        # with open(os.path.join(output_folder, re.sub(r'[><=]+', '_', prefix) + '.Release.' + f + '.pathway.txt'), 'w') as op2:
        #     op2.write('Pathway\tKept\n')
        #     for pwy in sorted(release_pathways_dat_pathways | release_pathways_dat_pathways_subset):
        #         if pwy in release_pathways_dat_pathways_subset:
        #             op2.write(pwy + '\tY\n' )
        #         else:
        #             op2.write(pwy + '\tN\n')

    # pgdb_helper_2(prefix + ' PreSAVI All_PreSAVI', rxn_list, all_pre_savi_reactions_dat_reaction_to_pathway,
    #               all_pre_savi_pathways_dat_reaction_to_pathway, all_pre_savi_pf_reactions)

    # pgdb_helper_2(prefix + ' Release All_Release', rxn_list, all_release_reactions_dat_reaction_to_pathway,
    #               all_release_pathways_dat_reaction_to_pathway, all_release_pf_reactions)


def main():
    prefix = 'Precision'
    over_size_rxn_with_ec_mapped_to_seq, over_size_jaccard_kept_rxn_to_seq, jaccard_kept_rxn_with_ec_mapped_to_seq \
        = filter_main(prefix)
    print()
    for f, pre_savi_reactions_dat_reaction_to_pathway, pre_savi_pathways_dat_reaction_to_pathway, \
        pre_savi_reactions_dat_reaction_to_pathway_subset, pre_savi_pathways_dat_reaction_to_pathway_subset, \
        release_reactions_dat_reaction_to_pathway, release_reactions_dat_reaction_to_pathway_subset, \
        release_pathways_dat_reaction_to_pathway_subset, release_pathways_dat_reaction_to_pathway in \
            pgdb_main('Size>=10', over_size_jaccard_kept_rxn_to_seq):
        print(f, len(pre_savi_reactions_dat_reaction_to_pathway), len(pre_savi_reactions_dat_reaction_to_pathway_subset),
              len(release_reactions_dat_reaction_to_pathway), len(release_pathways_dat_reaction_to_pathway_subset))
    print()
    for f, pre_savi_reactions_dat_reaction_to_pathway, pre_savi_pathways_dat_reaction_to_pathway, \
        pre_savi_reactions_dat_reaction_to_pathway_subset, pre_savi_pathways_dat_reaction_to_pathway_subset, \
        release_reactions_dat_reaction_to_pathway, release_reactions_dat_reaction_to_pathway_subset, \
        release_pathways_dat_reaction_to_pathway_subset, release_pathways_dat_reaction_to_pathway in \
        pgdb_main(prefix + '>=0.5', jaccard_kept_rxn_with_ec_mapped_to_seq):
        continue


if __name__ == '__main__':
    # main()
    prefix = 'Precision'
    _, over_size_jaccard_kept_rxn_to_seq, precision_jaccard_kept_rxn_with_ec_mapped_to_seq \
        = filter_main(prefix)

    prefix = 'Fscore'
    _, _, fscore_jaccard_kept_rxn_with_ec_mapped_to_seq \
        = filter_main(prefix)
    output_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21'
    with open(os.path.join(output_folder, 'size_10.rxn.txt'), 'w') as op1:
        for rxn in over_size_jaccard_kept_rxn_to_seq:
            op1.write(rxn + '\n')
    with open(os.path.join(output_folder, 'precision_0.5.rxn.txt'), 'w') as op2:
        for rxn in precision_jaccard_kept_rxn_with_ec_mapped_to_seq:
            op2.write(rxn + '\n')
    with open(os.path.join(output_folder, 'fscore_0.5.rxn.txt'), 'w') as op3:
        for rxn in fscore_jaccard_kept_rxn_with_ec_mapped_to_seq:
            op3.write(rxn + '\n')