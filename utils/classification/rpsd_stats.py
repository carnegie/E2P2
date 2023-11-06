import re

from Bio import SeqIO


def read_list_txt(txt_path, output_type='set'):
    if output_type == 'set':
        output = set()
    elif output_type == 'list':
        output = []
    else:
        return None
    with open(txt_path, 'r') as tp:
        for line in tp:
            line = line.rstrip()
            if output_type == 'set':
                output.add(line)
            elif output_type == 'list':
                output.append(line)
            else:
                continue
    return output


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


def filter_performances(performance_path, threshold=0.5):
    kept_classes = set()
    with open(performance_path, 'r') as pp:
        for line in pp:
            info = [i.strip() for i in line.split('\t')]
            try:
                input_class = info[0]
                fscore = float(info[-1])
                if fscore >= threshold:
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


def validate_rxn_to_ec(rxn_to_ec_dict):
    rxn_with_multiple_ec = set()
    ec_not_same_first_3 = {}
    ec_not_same_first_2 = {}
    ec_not_same_first_1 = {}
    for rxn in rxn_to_ec_dict:
        if len(rxn_to_ec_dict[rxn]) > 1:
            rxn_with_multiple_ec.add(rxn)
            first_3 = set(['.'.join(ec.split('.')[:3]) for ec in sorted(rxn_to_ec_dict[rxn])])
            if len(first_3) > 1:
                ec_not_same_first_3.setdefault(rxn, set(rxn_to_ec_dict[rxn]))
                first_2 = set(['.'.join(ec.split('.')[:2]) for ec in sorted(rxn_to_ec_dict[rxn])])
                if len(first_2) > 1:
                    ec_not_same_first_2.setdefault(rxn, set(rxn_to_ec_dict[rxn]))
                    first_1 = set(['.'.join(ec.split('.')[0]) for ec in sorted(rxn_to_ec_dict[rxn])])
                    if len(first_1) > 1:
                        ec_not_same_first_1.setdefault(rxn, set(rxn_to_ec_dict[rxn]))
    for key in ec_not_same_first_1:
        ec_not_same_first_2.pop(key, None)
        ec_not_same_first_3.pop(key, None)
    for key in ec_not_same_first_2:
        ec_not_same_first_3.pop(key, None)
    return rxn_with_multiple_ec, ec_not_same_first_1, ec_not_same_first_2, ec_not_same_first_3


def main():
    size = 10
    threshold = 0.5
    rpsd_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/release_2019-03-07/fasta/rpsd-4.2-20190307.fasta'
    ec_superseded_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/release_2019-03-07/maps/' \
                         'pf-EC-superseded.mapping'
    seq_to_annot, annot_to_seq = process_rpsd(rpsd_path, ec_superseded_path)
    print('RPSD Seq', len(seq_to_annot), 'RPSD Annot', len(annot_to_seq))
    ec_1_to_seq, ec_2_to_seq, ec_3_to_seq, ec_to_seq, rxn_to_seq = partition_annotations(annot_to_seq)
    print('RPSD EC', len(ec_to_seq), 'RPSD RXN', len(rxn_to_seq))

    official_ec_to_rxn_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/release_2019-03-07/maps/pf-official-EC-metacyc-RXN.mapping'
    rxn_to_ec_path = '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/release_2019-03-07/maps/pf-metacyc-RXN-EC.mapping'
    official_rxn_to_ec_dict, others_rxn_to_ec_dict = \
        read_metacyc_rxn_to_ec_maps(official_ec_to_rxn_path, rxn_to_ec_path)
    all_rxns = set(rxn_to_seq.keys())
    for ec in sorted(ec_to_seq.keys()):
        for rxn in official_rxn_to_ec_dict:
            if ec in official_rxn_to_ec_dict[rxn]:
                all_rxns.add(rxn)
        for rxn in others_rxn_to_ec_dict:
            if ec in others_rxn_to_ec_dict[rxn]:
                all_rxns.add(rxn)
    print('RPSD RXN w/ EC mapped', len(all_rxns))
    pgdb_plantcyc_reactions_path = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/reactions.plantcyc_02_16_21.ids.txt'
    pgdb_species_reactions_path = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/reactions.species_02_16_21.ids.txt'
    e2p2_reactions_path = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/e2p2.21.5-mapped.species_02_16_21.ids.txt'

    plantcyc_reactions = read_list_txt(pgdb_plantcyc_reactions_path)
    species_reactions = read_list_txt(pgdb_species_reactions_path)
    e2p2_reactions = read_list_txt(e2p2_reactions_path)

    print('Plantcyc RXN', len(plantcyc_reactions), 'RPSD Plantcyc RXN', len(all_rxns & plantcyc_reactions))
    print('Species RXN', len(species_reactions), 'RPSD Species RXN', len(all_rxns & species_reactions))
    print('E2P2 RXN', len(e2p2_reactions), 'RPSD E2P2 RXN', len(all_rxns & e2p2_reactions))
    rxn_seqs = set()
    for rxn in rxn_to_seq:
        rxn_seqs.update(rxn_to_seq[rxn])
    print('RPSD RXN Seqs', len(rxn_seqs))
    over_10_seqs = set()
    over_10_ec_set = set()
    for ec in ec_to_seq:
        if len(ec_to_seq[ec]) >= size:
            over_10_ec_set.add(ec)
            over_10_seqs.update(set(ec_to_seq[ec]))
    over_10_rxn_set = set()
    for rxn in rxn_to_seq:
        if len(rxn_to_seq[rxn]) >= size:
            over_10_rxn_set.add(rxn)
            over_10_seqs.update(set(rxn_to_seq[rxn]))

    over_10_all_rxn = set(over_10_rxn_set)
    for ec in sorted(over_10_ec_set):
        for rxn in official_rxn_to_ec_dict:
            if ec in official_rxn_to_ec_dict[rxn]:
                over_10_all_rxn.add(rxn)
        for rxn in others_rxn_to_ec_dict:
            if ec in others_rxn_to_ec_dict[rxn]:
                over_10_all_rxn.add(rxn)
    print('>=', size, 'Seq', len(over_10_seqs), '>=', size, 'RXN', len(over_10_rxn_set), '>=', size, 'EC', len(over_10_ec_set))
    print('>=', size, 'RPSD RXN w/ EC mapped', len(over_10_all_rxn))
    print('>=', size, 'Plantcyc RXN', len(over_10_all_rxn & plantcyc_reactions))
    print('>=', size, 'Species RXN', len(over_10_all_rxn & species_reactions))
    print('>=', size, 'E2P2 RXN', len(over_10_all_rxn & e2p2_reactions))
    jaccard_ef_mapping = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/clustering/merged/efclasses.jaccard.mapping'
    jaccard_ef_to_id, jacard_id_to_ef = read_ef_mapping(jaccard_ef_mapping, ec_superseded_path)

    jaccard_weight = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/clustering/merged/partitions/cross_validation/02_09_21/all/rpsd-4.2-20190307.ef.jaccard.blast.1e-2.label.weights.txt'
    jaccard_kept_classes = filter_performances(jaccard_weight, threshold=threshold)
    jaccard_kept_ids = set()
    for ef in sorted(jaccard_kept_classes):
        try:
            mapped_ids = jaccard_ef_to_id[ef]
            jaccard_kept_ids.update(mapped_ids)
            if len(mapped_ids) > 1:
                jaccard_ecs = set([i for i in sorted(mapped_ids) if re.search(r'^(\d+\.){3}\d+$', i)])
                jaccard_rxns = set([i for i in sorted(mapped_ids) if not re.search(r'^(\d+\.){3}\d+$', i)])
                # if len(jaccard_ecs) > 0 and len(jaccard_rxns) > 0:
                #     print(ef, sorted(jaccard_ecs), sorted(jaccard_rxns))
        except KeyError:
            print('EF not in map', ef)
    jaccard_kept_seqs = set()
    for annot in sorted(jaccard_kept_ids):
        try:
            seq_of_annot = annot_to_seq[annot]
            jaccard_kept_seqs.update(seq_of_annot)
        except KeyError:
            print('No seq found', annot)
    kept_jaccard_ecs = set([i for i in sorted(jaccard_kept_ids) if re.search(r'^(\d+\.){3}\d+$', i)])
    kept_jaccard_rxns = set([i for i in sorted(jaccard_kept_ids) if not re.search(r'^(\d+\.){3}\d+$', i)])
    print('Performance >=', threshold, 'Jaccard Class', len(jaccard_kept_classes), 'Jaccard IDs', len(jaccard_kept_ids),
          'Jaccard Seqs', len(jaccard_kept_seqs))
    print('Performance >=', threshold, 'Jaccard RXN', len(kept_jaccard_rxns), 'Jaccard EC', len(kept_jaccard_ecs))
    all_jaccard_rxns = set(kept_jaccard_rxns)
    for ec in sorted(kept_jaccard_ecs):
        for rxn in official_rxn_to_ec_dict:
            if ec in official_rxn_to_ec_dict[rxn]:
                all_jaccard_rxns.add(rxn)
        for rxn in others_rxn_to_ec_dict:
            if ec in others_rxn_to_ec_dict[rxn]:
                all_jaccard_rxns.add(rxn)
    print('Performance >=', threshold, 'Jaccard RXN w/ EC mapped', len(all_jaccard_rxns))
    print('Performance >=', threshold, 'Jaccard Plantcyc RXN', len(all_jaccard_rxns & plantcyc_reactions))
    print('Performance >=', threshold, 'Jaccard Species RXN', len(all_jaccard_rxns & species_reactions))
    print('Performance >=', threshold, 'Jaccard E2P2 RXN', len(all_jaccard_rxns & e2p2_reactions))

    official_rxn_to_ec_dict, others_rxn_to_ec_dict = \
        read_metacyc_rxn_to_ec_maps(official_ec_to_rxn_path, rxn_to_ec_path, ec_list=set(ec_to_seq.keys()),
                                    rxn_list=set(rxn_to_seq.keys()))
    print('RXN with no EC', len(set(rxn_to_seq.keys()) - set(official_rxn_to_ec_dict.keys()) - set(others_rxn_to_ec_dict.keys())))
    rxn_no_ec_seqs = set()
    for rxn in set(rxn_to_seq.keys()) - set(official_rxn_to_ec_dict.keys()) - set(others_rxn_to_ec_dict.keys()):
        rxn_no_ec_seqs.update(rxn_to_seq[rxn])
    print('RXN with no EC Seqs', len(rxn_no_ec_seqs))
    # for seq in sorted(rxn_no_ec_seqs):
    #     print(seq, seq_to_annot[seq])

    print('RXN with EC', len(set(official_rxn_to_ec_dict.keys()) | set(others_rxn_to_ec_dict.keys())))
    rxn_with_ec_seqs = set()
    for rxn in set(official_rxn_to_ec_dict.keys()) | set(others_rxn_to_ec_dict.keys()):
        rxn_with_ec_seqs.update(rxn_to_seq[rxn])
    print('RXN with EC Seqs', len(rxn_with_ec_seqs))

    print('Official RXN with EC', len(official_rxn_to_ec_dict))
    rxn_with_official_ec_seqs = set()
    for rxn in official_rxn_to_ec_dict:
        rxn_with_official_ec_seqs.update(rxn_to_seq[rxn])
    print('Official RXN with EC Seqs', len(rxn_with_official_ec_seqs))

    print('Non official RXN with EC', len(others_rxn_to_ec_dict))
    rxn_with_no_official_ec_seqs = set()
    for rxn in others_rxn_to_ec_dict:
        rxn_with_no_official_ec_seqs.update(rxn_to_seq[rxn])
    print('Non Official RXN with EC Seqs', len(rxn_with_no_official_ec_seqs))

    official_rxn_with_multiple_ec, official_ec_not_same_first_1, official_ec_not_same_first_2, official_ec_not_same_first_3 = \
        validate_rxn_to_ec(official_rxn_to_ec_dict)
    print('Official RXN with multiple EC', len(official_rxn_with_multiple_ec))
    official_rxn_with_multiple_ec_seqs = set()
    for rxn in official_rxn_with_multiple_ec:
        official_rxn_with_multiple_ec_seqs.update(rxn_to_seq[rxn])
    print('Official RXN with multiple EC Seqs', len(official_rxn_with_multiple_ec_seqs))

    print('Official RXN with multiple EC, third part not the same', len(official_ec_not_same_first_3))
    print('Official RXN with multiple EC, second part not the same', len(official_ec_not_same_first_2))
    print('Official RXN with multiple EC, first part not the same', len(official_ec_not_same_first_1))
    non_rxn_with_multiple_ec, non_ec_not_same_first_1, non_ec_not_same_first_2, non_ec_not_same_first_3 = \
        validate_rxn_to_ec(others_rxn_to_ec_dict)
    print('Non Official RXN with multiple EC', len(non_rxn_with_multiple_ec))
    non_rxn_with_multiple_ec_seqs = set()
    for rxn in non_rxn_with_multiple_ec:
        non_rxn_with_multiple_ec_seqs.update(rxn_to_seq[rxn])
    print('Non Official RXN with multiple EC Seqs', len(non_rxn_with_multiple_ec_seqs))

    print('Non Official RXN with multiple EC, third part not the same', len(non_ec_not_same_first_3))
    print('Non Official RXN with multiple EC, second part not the same', len(non_ec_not_same_first_2))
    print('Non Official RXN with multiple EC, first part not the same', len(non_ec_not_same_first_1))

    rxn_in_official_and_others = {}
    for rxn in set(others_rxn_to_ec_dict.keys()) & set(official_rxn_to_ec_dict.keys()):
        rxn_in_official_and_others.setdefault(rxn, set(others_rxn_to_ec_dict[rxn]) | set(official_rxn_to_ec_dict[rxn]))
    both_rxn_with_multiple_ec, both_ec_not_same_first_1, both_ec_not_same_first_2, both_ec_not_same_first_3 = \
        validate_rxn_to_ec(rxn_in_official_and_others)
    print('Both RXN with EC', len(rxn_in_official_and_others))
    print('Both RXN with multiple EC', len(both_rxn_with_multiple_ec))

    both_rxn_with_multiple_ec_seqs = set()
    for rxn in both_rxn_with_multiple_ec:
        both_rxn_with_multiple_ec_seqs.update(rxn_to_seq[rxn])
    print('Both RXN with multiple EC Seqs', len(both_rxn_with_multiple_ec_seqs))

    print('Both RXN with multiple EC, third part not the same', len(both_ec_not_same_first_3))
    print('Both RXN with multiple EC, second part not the same', len(both_ec_not_same_first_2))
    print('Both RXN with multiple EC, first part not the same', len(both_ec_not_same_first_1))



    # print(plantcyc_reactions)