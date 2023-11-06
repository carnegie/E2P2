import os

from utils.classification.pgdb_filter_postprocess import reverse_dict
from utils.classification.pgdb_filter_test import process_pgdb_folder

pgdb_list = ["aracyc", "brapa_fpsccyc", "chlamycyc", "corncyc", "mtruncatulacyc", "oryzacyc", "poplarcyc",
             "quinoacyc", "soycyc", "tomatocyc"]


def get_pathways_from_savi_txt(savi_text_path):
    pathways = set()
    with open(savi_text_path, 'r') as fp:
        for line in fp:
            if line.startswith('PWY-ID\t'):
                continue
            info = [i.strip() for i in line.split('\t')]
            pathways.add(info[0])
    return pathways


def process_savi_folder(savi_folder):
    accepted_txt = os.path.join(savi_folder, 'Accepted.txt')
    manual_validate_txt = os.path.join(savi_folder, 'Manual-to-validate.txt')
    if os.path.isfile(accepted_txt) and os.path.isfile(manual_validate_txt):
        accepted_pathways = get_pathways_from_savi_txt(accepted_txt)
        manual_validate_pathways = get_pathways_from_savi_txt(manual_validate_txt)
        return accepted_pathways | manual_validate_pathways
    else:
        print('SAVI folder process error', savi_folder)
        raise SystemError


def write_removed_pathways(release_pathways, kept_pathways, output_path):
    with open(output_path, 'w') as op:
        for pwy in sorted(set(release_pathways) - set(kept_pathways)):
            op.write(pwy + '\n')


if __name__ == '__main__':
    release_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_release_02_16_21'
    savi_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/savi_filtered_02_16_21'
    output_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/savi_removed_02_16_21'
    for pgdb in pgdb_list:
        release_pgdb_folder = os.path.join(release_folder, pgdb)
        _, release_pathways_dat_reaction_to_pathway, _ = process_pgdb_folder(release_pgdb_folder)
        release_pathways_dat_pathway_to_reaction = reverse_dict(release_pathways_dat_reaction_to_pathway)

        fscore_folder_name = pgdb + '_fscore'
        precision_folder_name = pgdb + '_precision'
        size_folder_name = pgdb + '_size'

        fscore_savi_path = os.path.join(savi_folder, fscore_folder_name)
        precision_savi_path = os.path.join(savi_folder, precision_folder_name)
        size_savi_path = os.path.join(savi_folder, size_folder_name)

        fscore_pathways = process_savi_folder(fscore_savi_path)
        precision_pathways = process_savi_folder(precision_savi_path)
        size_pathways = process_savi_folder(size_savi_path)

        fscore_output = os.path.join(output_folder, fscore_folder_name + '.txt')
        precision_output = os.path.join(output_folder, precision_folder_name + '.txt')
        size_output = os.path.join(output_folder, size_folder_name + '.txt')

        write_removed_pathways(set(release_pathways_dat_pathway_to_reaction.keys()), fscore_pathways, fscore_output)
        write_removed_pathways(set(release_pathways_dat_pathway_to_reaction.keys()), precision_pathways, precision_output)
        write_removed_pathways(set(release_pathways_dat_pathway_to_reaction.keys()), size_pathways, size_output)