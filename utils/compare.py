import re


def read_priam_sequence_ec_itr(sequence_ecs_fp, split=r'\s+|\|', skip='#'):
    query_id, priam_results = '', []
    for line in sequence_ecs_fp:
        line = line.strip()
        if line.startswith('>'):
            if len(query_id) > 0:
                yield query_id, priam_results
            query_id = re.sub(r'^>', '', re.split(split, line)[0])
            query_id, priam_results = query_id, []
        else:
            if not line.startswith(skip) and len(line) > 0:
                try:
                    info = line.split('\t')
                    ef_class = info[0].strip()
                    try:
                        e_value = float(info[2])
                        priam_results.append((ef_class, e_value))
                    except ValueError:
                        continue
                except IndexError:
                    continue
            else:
                continue
    if len(query_id) > 0:
        yield query_id, priam_results


def read_map_file(file_path, key_idx=0, val_idx=1):
    mapping_dict = {}
    with open(file_path, 'r') as fp:
        for line in fp:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            info = line.strip().split('\t')
            try:
                mapping_dict[info[key_idx]].append(info[val_idx])
            except KeyError:
                mapping_dict.setdefault(info[key_idx], [info[val_idx]])
            except IndexError:
                continue
    return mapping_dict


def read_priam_output(path_to_sequence_ecs, split=r'\s+|\|', skip='#', efmap=None):
    priam_output_dict = {}
    with open(path_to_sequence_ecs) as ptse:
        for query_id, priam_results in read_priam_sequence_ec_itr(ptse, split=split, skip=skip):
            if efmap:
                mapped_priam_results = []
                for res in priam_results:
                    res_cls = res[0]
                    res_score = res[1]
                    try:
                        mapped_cls = efmap[res_cls.replace('#', '')]
                        if res_cls.startswith('#'):
                            mapped_cls = ['#' + cls for cls in mapped_cls]
                        for cls in mapped_cls:
                            mapped_priam_results.append((cls, res_score))
                    except KeyError:
                        print("Not mapped", res_cls)
                priam_results = mapped_priam_results

            priam_output_dict.setdefault(query_id, priam_results)
    return priam_output_dict


def convert_list_of_tuples_to_dict(list_of_tuples, key_idx=0, val_idx=1):
    output_dict = {}
    for tup in list_of_tuples:
        try:
            key = tup[key_idx]
            val = tup[val_idx]
            try:
                output_dict[key].add(val)
            except KeyError:
                output_dict.setdefault(key, {val})
        except IndexError:
            continue
    return output_dict


def compare_og_to_new(priam_f0_dict, priam_test_dict):
    for key in sorted(set(priam_f0_dict.keys()) - set(priam_test_dict.keys())):
        print('In OG only', key, priam_f0_dict[key])

    for key in sorted(set(priam_test_dict.keys()) - set(priam_f0_dict.keys())):
        print('In 5 CPU only', key, priam_test_dict[key])

    for key in sorted(set(priam_test_dict.keys()) & set(priam_f0_dict.keys())):
        res_test = convert_list_of_tuples_to_dict(priam_test_dict[key])
        res_f0 = convert_list_of_tuples_to_dict(priam_f0_dict[key])
        same = True
        if len(res_test) != len(res_f0):
            same = False
        elif set(res_test.keys()) != set(res_f0.keys()):
            same = False
        else:
            for res_key in set(res_test.keys()) & set(res_f0.keys()):
                if res_test[res_key] != res_f0[res_key]:
                    same = False
        if same is False:
            print('Diff Res, 5 CPU', priam_test_dict[key], 'OG', priam_f0_dict[key])


if __name__ == '__main__':
    priam_f0_path = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/priam/cpu_compare/Priam-f0/' \
                    'sequenceECs.txt'
    priam_test_path = '/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-e2p2/priam/cpu_compare/Priam-test/' \
                      'sequenceECs.txt'
    ef_class_map_path = '/Users/bxue/Documents/Carnegie/SourceData/PMN/RPSD/release_2019-03-07/maps/efclasses.mapping'
    ef_class_map = read_map_file(ef_class_map_path)
    print("All Res")
    priam_f0_dict = read_priam_output(priam_f0_path, split='#_#', efmap=ef_class_map, skip='$')
    priam_test_dict = read_priam_output(priam_test_path, split='#_#', efmap=ef_class_map, skip='$')
    compare_og_to_new(priam_f0_dict, priam_test_dict)

    print("Skiped '#'")
    priam_f0_dict = read_priam_output(priam_f0_path, split='#_#', efmap=ef_class_map, skip='#')
    priam_test_dict = read_priam_output(priam_test_path, split='#_#', efmap=ef_class_map, skip='#')
    compare_og_to_new(priam_f0_dict, priam_test_dict)
