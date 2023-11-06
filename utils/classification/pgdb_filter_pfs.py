import os
import re
from shutil import copyfile

from utils.classification.pgdb_filter_postprocess import read_rxn_list
from utils.classification.pgdb_filter_test import read_metacyc_rxn_to_ec_maps


def read_pf_itr(fp):
    seq_id, attrs = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith("ID"):
            if seq_id:
                yield seq_id, attrs
            seq_id, attrs = re.sub(r'^ID\t', '', line), []
        else:
            info = [i.strip() for i in re.split(r'\t', line)]
            try:
                attrs.append((info[0], info[1]))
            except IndexError:
                continue
    if seq_id:
        yield seq_id, attrs


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


def filter_pf_file(pf_path, rxn_list, output_prefix):
    output_path = os.path.splitext(pf_path)[0] + '.' + output_prefix + os.path.splitext(pf_path)[1]
    bak_pf_path = pf_path + '.bak'
    copyfile(pf_path, bak_pf_path)
    with open(bak_pf_path, 'r') as bpp, open(output_path, 'w') as pp:
        for seq_id, attrs in read_pf_itr(bpp):
            attrs_of_seq_id = []
            metacycs_of_seq_id = []
            for attr in attrs:
                try:
                    attr_name = attr[0]
                    attr_val = attr[1]
                    if attr_name == 'METACYC' and attr_val in rxn_list:
                        metacycs_of_seq_id.append(attr)
                    elif attr_name != 'METACYC':
                        attrs_of_seq_id.append(attr)
                except IndexError:
                    print('Attr Error', attr)
            if len(metacycs_of_seq_id) == 0:
                continue
            else:
                pp.write('ID\t' + seq_id + '\n')
                for attr in attrs_of_seq_id + metacycs_of_seq_id:
                    pp.write('\t'.join(attr) + '\n')
                pp.write('//\n')


if __name__ == '__main__':
    # brapa_fpsccyc, mtruncatulacyc, poplarcyc, qunioacyc, tomatocyc
    pf_path = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_input_02_16_21/tomatocyc/ITAG3.2_proteins.fa.e2p2v4.orxn.revised.pf'

    fscore_rxn_list = read_rxn_list('/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/rxn.fscore_0.5.txt')
    precision_rxn_list = read_rxn_list(
        '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/rxn.precision_0.5.txt')
    size_rxn_list = read_rxn_list('/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/rxn.size_10.txt')

    filter_pf_file(pf_path, fscore_rxn_list, 'fscore')
    filter_pf_file(pf_path, precision_rxn_list, 'precision')
    filter_pf_file(pf_path, size_rxn_list, 'size')
