import os


def read_ef_map(file_path):
    id_to_ef_dict = {}
    ef_to_id_dict = {}
    with open(file_path, 'r') as fp:
        for line in fp:
            info = line.rstrip().split('\t')
            id_to_ef_dict.setdefault(info[1], info[0])
            ef_to_id_dict.setdefault(info[0], info[1])
    return id_to_ef_dict, ef_to_id_dict


def add_value_list_to_key_dict(value_list, key, input_dict):
    try:
        input_dict[key].update(value_list)
    except KeyError:
        input_dict.setdefault(key, set(value_list))


def read_seta_annotation(file_path):
    ec_no_rxn_set = set()
    rxn_no_ec_set = set()
    rxn_to_ec_dict = {}
    ec_to_rxn_dict = {}
    prot_to_ef_dict = {}
    ef_to_prot_dict = {}
    with open(file_path, 'r') as fp:
        for line in fp:
            line = line.rstrip()
            info = line.split('|')
            uniprot_acc = info[0]
            function_ids = info[1].split(',')
            for f_id in function_ids:
                if f_id.startswith("NA:"):
                    ef_id = f_id.replace('NA:', '')
                    ec_no_rxn_set.add(ef_id)
                    add_value_list_to_key_dict([ef_id], uniprot_acc, prot_to_ef_dict)
                    add_value_list_to_key_dict([uniprot_acc], ef_id, ef_to_prot_dict)
                elif f_id.endswith(':NA'):
                    ef_id = f_id.replace(':NA', '')
                    rxn_no_ec_set.add(ef_id)
                    add_value_list_to_key_dict([ef_id], uniprot_acc, prot_to_ef_dict)
                    add_value_list_to_key_dict([uniprot_acc], ef_id, ef_to_prot_dict)
                else:
                    rxn_ec_pair = [f.strip() for f in f_id.split(':') if f != '']
                    if len(rxn_ec_pair) >= 2:
                        add_value_list_to_key_dict(rxn_ec_pair, uniprot_acc, prot_to_ef_dict)
                        add_value_list_to_key_dict([rxn_ec_pair[1]], rxn_ec_pair[0], rxn_to_ec_dict)
                        add_value_list_to_key_dict([rxn_ec_pair[0]], rxn_ec_pair[1], ec_to_rxn_dict)
                        add_value_list_to_key_dict([uniprot_acc], rxn_ec_pair[0], ef_to_prot_dict)
                        add_value_list_to_key_dict([uniprot_acc], rxn_ec_pair[1], ef_to_prot_dict)
                    else:
                        print("Error Pair", rxn_ec_pair)
    ec_with_rxn_and_no_rxn_set = ec_no_rxn_set & set(ec_to_rxn_dict.keys())
    ec_no_rxn_set -= set(ec_to_rxn_dict.keys())
    ec_no_rxn_dict = {k:set() for k in ec_no_rxn_set}
    rxn_with_ec_and_no_ec_set = rxn_no_ec_set & set(rxn_to_ec_dict.keys())
    rxn_no_ec_set -= set(rxn_to_ec_dict.keys())
    rxn_no_ec_dict = {k: set() for k in rxn_no_ec_set}
    return ec_no_rxn_dict, rxn_no_ec_dict, rxn_to_ec_dict, ec_to_rxn_dict, prot_to_ef_dict, ef_to_prot_dict, \
           ec_with_rxn_and_no_rxn_set, rxn_with_ec_and_no_ec_set


def partition_dict_by_value_size(map_dict_1, map_1_keys_to_plus_one, map_dict_2, map_2_keys_to_plus_one):
    map_1_to_1_dict = {}
    map_1_to_1_to_n_dict = {}
    map_1_to_n_dict = {}
    for key_1 in sorted(map_dict_1.keys()):
        values_of_key_1 = map_dict_1[key_1]
        len_values_of_key_1 = len(values_of_key_1)
        if key_1 in map_1_keys_to_plus_one:
            len_values_of_key_1 += 1
        if len_values_of_key_1 == 1:
            value_of_key = sorted(values_of_key_1)[0]
            try:
                values_of_value_of_key = map_dict_2[value_of_key]
                len_values_of_value_of_key = len(values_of_value_of_key)
                if value_of_key in map_2_keys_to_plus_one:
                    len_values_of_value_of_key += 1
                if len_values_of_value_of_key == 1:
                    map_1_to_1_dict.setdefault(key_1, values_of_key_1)
                else:
                    map_1_to_1_to_n_dict.setdefault(key_1, values_of_key_1)
            except KeyError:
                print("Values of", value_of_key, "not found.")
        else:
            map_1_to_n_dict.setdefault(key_1, values_of_key_1)
    return map_1_to_1_dict, map_1_to_1_to_n_dict, map_1_to_n_dict


def group_rxn_ec(rxn_to_ec_dict, ec_to_rxn_dict, ec_with_rxn_and_no_rxn_set, rxn_with_ec_and_no_ec_set):
    rxn_1_to_1_ec, rxn_1_to_1_ec_to_n_rxn, rxn_1_to_n_ec = partition_dict_by_value_size(rxn_to_ec_dict,
                                                                                        rxn_with_ec_and_no_ec_set,
                                                                                        ec_to_rxn_dict,
                                                                                        ec_with_rxn_and_no_rxn_set)
    ec_1_to_1_rxn, ec_1_to_1_rxn_to_n_ec, ec_1_to_n_rxn = partition_dict_by_value_size(ec_to_rxn_dict,
                                                                                       ec_with_rxn_and_no_rxn_set,
                                                                                       rxn_to_ec_dict,
                                                                                       rxn_with_ec_and_no_ec_set)
    return rxn_1_to_1_ec, rxn_1_to_1_ec_to_n_rxn, rxn_1_to_n_ec, ec_1_to_1_rxn, ec_1_to_1_rxn_to_n_ec, ec_1_to_n_rxn


def read_performance_output(file_path, efs_to_ids):
    label_performance = {}
    instance_performance = {}
    write_label_flag = False
    write_instance_flag = False
    with open(file_path, 'r') as fp:
        for line in fp:
            if line.startswith('#Performance of each label'):
                write_label_flag = True
            if line.startswith('#Performance of each instance'):
                write_label_flag = False
                write_instance_flag = True
            if write_label_flag and not line.startswith('#'):
                try:
                    info = line.rstrip().split('\t')
                    ef_class = info[0]
                    precision = info[-3]
                    recall = info[-2]
                    f_score = info[-1]
                    try:
                        mapped_id = efs_to_ids[ef_class]
                        label_performance.setdefault(mapped_id, (precision, recall, f_score))
                    except KeyError:
                        print('Error EF', ef_class, info)
                        continue
                except IndexError:
                    print("Index Error", info)
                    continue
            if write_instance_flag and not line.startswith('#'):
                try:
                    info = line.rstrip().split('\t')
                    uniprot_id = info[0]
                    precision = info[-3]
                    recall = info[-2]
                    f_score = info[-1]
                    instance_performance.setdefault(uniprot_id, (precision, recall, f_score))
                except IndexError:
                    print("Index Error", info)
                    continue
    return label_performance, instance_performance


def get_ef_info(func_id, id_to_ef_dict, id_to_prot_dict, label_performance):
    try:
        ef_of_id = id_to_ef_dict[func_id]
    except KeyError:
        ef_of_id = "NA"
        print('ID', func_id, 'not found in ID to EF map')
    try:
        prots_of_id = id_to_prot_dict[func_id]
    except KeyError:
        prots_of_id = set()
        print('ID', func_id, 'not found in ID to Protein map')
    try:
        performances_of_id = label_performance[func_id]
    except KeyError:
        performances_of_id = ("NA", "NA", "NA")
        # print('ID', func_id, 'not found in ID to Performance map')
    return ef_of_id, prots_of_id, performances_of_id


def func_id_info_iterator(id_map_dict, id_to_ef_dict, id_to_prot_dict, label_performance, instance_performance):
    for func_id in sorted(id_map_dict.keys()):
        mapped_ids = id_map_dict[func_id]
        ef_id, prots_of_id, performances_of_id = get_ef_info(func_id, id_to_ef_dict, id_to_prot_dict, label_performance)
        prots_of_id_in_test = set(prots_of_id) & set(instance_performance.keys())
        mapped_id_info = []
        for mapped_func_id in sorted(mapped_ids):
            mapped_ef_id, prots_of_mapped_id, performances_of_mapped_id = get_ef_info(mapped_func_id, id_to_ef_dict,
                                                                                      id_to_prot_dict, label_performance)
            prots_of_mapped_id_in_test = set(prots_of_mapped_id) & set(instance_performance.keys())
            mapped_id_info.append((mapped_ef_id, prots_of_mapped_id, prots_of_mapped_id_in_test,
                                   performances_of_mapped_id))
        yield func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info


def calculate_average_performance_tuples(combined_performances):
    na_removed = [t for t in combined_performances if 'NA' not in t]
    if len(na_removed) == 0:
        return "NA", "NA", "NA"
    else:
        num_of_performance = len(na_removed)
        sum_precision = sum([float(t[0]) for t in na_removed])
        sum_recall = sum([float(t[1]) for t in na_removed])
        sum_f_score = sum([float(t[2]) for t in na_removed])
        return str(sum_precision/num_of_performance), str(sum_recall/num_of_performance), str(sum_f_score/num_of_performance)


def test():
    seta_annotation_path = '../test/setA-annotation'
    ef_map_path = '../test/setA-efmaps'
    performance_path = '../test/rpsd-4.2-20190307.e2p2v4.performance.txt'
    output_folder = '../test/performances'

    ec_no_rxn_dict, rxn_no_ec_dict, rxn_to_ec_dict, ec_to_rxn_dict, prot_to_ef_dict, ef_to_prot_dict, \
        ec_with_rxn_and_no_rxn_set, rxn_with_ec_and_no_ec_set = read_seta_annotation(seta_annotation_path)

    id_to_ef_dict, ef_to_id_dict = read_ef_map(ef_map_path)
    label_performance, instance_performance = read_performance_output(performance_path, ef_to_id_dict)
    rxn_1_to_1_ec, rxn_1_to_1_ec_to_n_rxn, rxn_1_to_n_ec, ec_1_to_1_rxn, ec_1_to_1_rxn_to_n_ec, ec_1_to_n_rxn = \
        group_rxn_ec(rxn_to_ec_dict, ec_to_rxn_dict, ec_with_rxn_and_no_rxn_set, rxn_with_ec_and_no_ec_set)

    original_max_ef_prot_to_uniprot = {}
    combined_max_ef_prot_to_uniprot = {}

    for uniprot_id in sorted(prot_to_ef_dict.keys()):
        ef_classes_of_uniprot = {k:ef_to_prot_dict[k] for k in prot_to_ef_dict[uniprot_id]}
        combined_ef_classes_of_uniprot = {}
        for ef_class in sorted(prot_to_ef_dict[uniprot_id]):
            try:
                mapped_ef_classes = ec_1_to_1_rxn[ef_class]
                for mapped_ef in sorted(mapped_ef_classes):
                    add_value_list_to_key_dict(ef_to_prot_dict[mapped_ef], mapped_ef,
                                                combined_ef_classes_of_uniprot)
            except KeyError:
                add_value_list_to_key_dict(ef_to_prot_dict[ef_class], ef_class,
                                            combined_ef_classes_of_uniprot)
        original_ef_prot_lens = []
        combined_ef_prot_lens = []
        for ef_class in sorted(ef_classes_of_uniprot.keys()):
            original_ef_prot_lens.append(len(ef_classes_of_uniprot[ef_class]))
        for ef_class in sorted(combined_ef_classes_of_uniprot.keys()):
            combined_ef_prot_lens.append(len(combined_ef_classes_of_uniprot[ef_class]))

        try:
            original_max_ef_prot_to_uniprot[max(original_ef_prot_lens)].add(uniprot_id)
        except KeyError:
            original_max_ef_prot_to_uniprot.setdefault(max(original_ef_prot_lens), {uniprot_id})
        try:
            combined_max_ef_prot_to_uniprot[max(combined_ef_prot_lens)].add(uniprot_id)
        except KeyError:
            combined_max_ef_prot_to_uniprot.setdefault(max(combined_ef_prot_lens), {uniprot_id})
    max_ef_prots_len_to_num_of_uniprots_path = os.path.join(output_folder, "EFclasses_to_prots_max_to_prots_len.tsv")
    with open(max_ef_prots_len_to_num_of_uniprots_path, 'w') as op:
        op.write('Max_Num_of_prots_of_EFs\tNum_of_prots\tPercentage\n')
        for key in sorted(original_max_ef_prot_to_uniprot.keys()):
            op.write(str(key) + '\t' + str(len(original_max_ef_prot_to_uniprot[key])) + '\t' +
                     str(len(original_max_ef_prot_to_uniprot[key])/57514) + '\n')

    print('RXN 1 to 1 EC', len(rxn_1_to_1_ec))
    print('RXN 1 to 1 EC to N RXN', len(rxn_1_to_1_ec_to_n_rxn))
    print('RXN 1 to N EC', len(rxn_1_to_n_ec))
    print('RXN 1 to 0 EC', len(rxn_no_ec_dict))
    print('EC 1 to 1 RXN', len(ec_1_to_1_rxn))
    print('EC 1 to 1 RXN to N EC', len(ec_1_to_1_rxn_to_n_ec))
    print('EC 1 to N RXN', len(ec_1_to_n_rxn))
    print('EC 1 to 0 RXN', len(ec_no_rxn_dict))
    print('Total', len(rxn_1_to_1_ec) + len(rxn_1_to_1_ec_to_n_rxn) + len(rxn_1_to_n_ec) + len(rxn_no_ec_dict) +
          len(ec_1_to_1_rxn) + len(ec_1_to_1_rxn_to_n_ec) + len(ec_1_to_n_rxn) + len(ec_no_rxn_dict))
    print('Total (combine 1 to 1)', len(rxn_1_to_1_ec) + len(rxn_1_to_1_ec_to_n_rxn) + len(rxn_1_to_n_ec) +
          len(rxn_no_ec_dict) + len(ec_1_to_1_rxn_to_n_ec) + len(ec_1_to_n_rxn) + len(ec_no_rxn_dict))
    combine_1_to_1_path = os.path.join(output_folder, 'EFclasses_performances.combine_1_to_1.tsv')
    original_path = os.path.join(output_folder, 'EFclasses_performances.original.tsv')
    with open(combine_1_to_1_path, 'w') as fp_1_to_1, open(original_path, 'w') as op:
        fp_1_to_1.write('ID\tEFClass\tNumOfSeq(RPSD)\tNumOfSeq(Test)\tPrecision\tRecall\tFscore\tType\tGroup\n')
        op.write('ID\tEFClass\tNumOfSeq(RPSD)\tNumOfSeq(Test)\tPrecision\tRecall\tFscore\tType\tGroup\n')
        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(rxn_1_to_1_ec, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "RXN", "RXN_1_to_1_EC"]) + '\n')

            combined_prots_of_id = set()
            combined_prots_of_id.update(prots_of_id)
            combined_prots_of_id_in_test = set()
            combined_prots_of_id_in_test.update(prots_of_id_in_test)
            combined_performances_of_id = [performances_of_id]
            for mapped_tuple in mapped_id_info:
                prots_of_mapped_id = mapped_tuple[1]
                prots_of_mapped_id_in_test = mapped_tuple[2]
                combined_prots_of_id.update(prots_of_mapped_id)
                combined_prots_of_id_in_test.update(prots_of_mapped_id_in_test)
                combined_performances_of_id.append(mapped_tuple[3])
            average_performances_of_id = calculate_average_performance_tuples(combined_performances_of_id)
            fp_1_to_1.write('\t'.join([func_id, ef_id, str(len(combined_prots_of_id)),
                                       str(len(combined_prots_of_id_in_test)),
                                       '\t'.join(average_performances_of_id), "RXN", "RXN_1_to_1_EC"]) + '\n')
        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(rxn_1_to_1_ec_to_n_rxn, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            fp_1_to_1.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                       '\t'.join(performances_of_id), "RXN", "RXN_1_to_1_EC_to_N_RXN"]) + '\n')
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "RXN", "RXN_1_to_1_EC_to_N_RXN"]) + '\n')
        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(rxn_1_to_n_ec, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            fp_1_to_1.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                       '\t'.join(performances_of_id), "RXN", "RXN_1_to_N_EC"]) + '\n')
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "RXN", "RXN_1_to_N_EC"]) + '\n')
        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(rxn_no_ec_dict, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            fp_1_to_1.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                       '\t'.join(performances_of_id), "RXN", "RXN_1_to_0_EC"]) + '\n')
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "RXN", "RXN_1_to_0_EC"]) + '\n')

        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(ec_1_to_1_rxn, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "EC", "EC_1_to_1_RXN"]) + '\n')
        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(ec_1_to_1_rxn_to_n_ec, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            fp_1_to_1.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                       '\t'.join(performances_of_id), "EC", "EC_1_to_1_RXN_to_N_EC"]) + '\n')
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "EC", "EC_1_to_1_RXN_to_N_EC"]) + '\n')
        with open('../test/performances/EC_1_to_N_RXN.tsv', 'w') as tsvop:
            for ec in sorted(ec_1_to_n_rxn.keys()):
                if len(ec_1_to_n_rxn[ec]) > 1:
                    tsvop.write(ec + '\t' + ','.join(sorted(ec_1_to_n_rxn[ec])) + '\n')
        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(ec_1_to_n_rxn, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            fp_1_to_1.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                       '\t'.join(performances_of_id), "EC", "EC_1_to_N_RXN"]) + '\n')
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "EC", "EC_1_to_N_RXN"]) + '\n')
        for func_id, ef_id, prots_of_id, prots_of_id_in_test, performances_of_id, mapped_id_info in \
                func_id_info_iterator(ec_no_rxn_dict, id_to_ef_dict, ef_to_prot_dict, label_performance,
                                      instance_performance):
            fp_1_to_1.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                       '\t'.join(performances_of_id), "EC", "EC_1_to_0_RXN"]) + '\n')
            op.write('\t'.join([func_id, ef_id, str(len(prots_of_id)), str(len(prots_of_id_in_test)),
                                '\t'.join(performances_of_id), "EC", "EC_1_to_0_RXN"]) + '\n')


if __name__ == '__main__':
    test()
