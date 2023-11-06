import os

from utils.classification.pgdb_filter_postprocess import reverse_dict
from utils.classification.pgdb_filter_test import process_pgdb_folder, process_reactions_dat, process_pathways_dat, \
    process_pf_file


def process_pathologic_folder(pgdb_folder_path):
    reactions_dat_reaction_to_pathway = {}
    pathways_dat_reaction_to_pathway = {}
    pf_reactions = set()
    for f in os.listdir(pgdb_folder_path):
        file_path = os.path.join(pgdb_folder_path, f)
        if os.path.isdir(file_path) and f == 'data':
            reactions_dat_path = os.path.join(file_path, 'reactions.dat')
            pathways_dat_path = os.path.join(file_path, 'pathways.dat')
            if os.path.isfile(reactions_dat_path) and os.path.isfile(pathways_dat_path):
                reactions_dat_reaction_to_pathway = process_reactions_dat(reactions_dat_path)
                pathways_dat_reaction_to_pathway = process_pathways_dat(pathways_dat_path)
            else:
                print('Data Folder error', file_path)
                raise SystemError
        elif os.path.isfile(file_path) and file_path.endswith('pf'):
            pf_reactions = process_pf_file(file_path)
            if len(pf_reactions) == 0:
                print('pf file error', file_path)
                raise SystemError
        else:
            continue
    return reactions_dat_reaction_to_pathway, pathways_dat_reaction_to_pathway, pf_reactions


def process_main(pgdb, pre_savi_folder, pathologic_folder, release_folder):
    pre_savi_pgdb_folder = os.path.join(pre_savi_folder, pgdb)
    release_pgdb_folder = os.path.join(release_folder, pgdb)

    fscore_patho_folder = os.path.join(pathologic_folder, 'fscore', pgdb)
    precision_patho_folder = os.path.join(pathologic_folder, 'precision', pgdb)
    size_patho_folder = os.path.join(pathologic_folder, 'size', pgdb)


    _, pre_savi_pathways_dat_reaction_to_pathway, pre_savi_pf_reactions = process_pgdb_folder(pre_savi_pgdb_folder)
    _, release_pathways_dat_reaction_to_pathway, release_pf_reactions = process_pgdb_folder(release_pgdb_folder)

    _, fscore_pathways_dat_reaction_to_pathway, fscore_pf_reactions = process_pathologic_folder(fscore_patho_folder)
    _, precision_pathways_dat_reaction_to_pathway, precision_pf_reactions = process_pathologic_folder(precision_patho_folder)
    _, size_pathways_dat_reaction_to_pathway, size_pf_reactions = process_pathologic_folder(size_patho_folder)

    return \
        pre_savi_pathways_dat_reaction_to_pathway, pre_savi_pf_reactions, \
        release_pathways_dat_reaction_to_pathway, release_pf_reactions, \
        fscore_pathways_dat_reaction_to_pathway, fscore_pf_reactions, \
        precision_pathways_dat_reaction_to_pathway, precision_pf_reactions, \
        size_pathways_dat_reaction_to_pathway, size_pf_reactions


if __name__ == '__main__':
    # pgdb_list = ['aracyc', 'chlamycyc', 'corncyc', 'oryzacyc', 'soycyc']
    pgdb_list = ["aracyc", "brapa_fpsccyc", "chlamycyc", "corncyc", "mtruncatulacyc", "oryzacyc", "poplarcyc",
                 "quinoacyc", "soycyc", "tomatocyc"]

    pre_savi_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_preSAVI_02_16_21'
    release_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_release_02_16_21'
    pathologic_folder = '/Users/bxue/Documents/Carnegie/PMNProject/pmn_data_02_16_21/pgdbs_filtered_02_16_21'

    for pgdb in pgdb_list:
        pre_savi_pathways_dat_reaction_to_pathway, pre_savi_pf_reactions, \
        release_pathways_dat_reaction_to_pathway, release_pf_reactions, \
        fscore_pathways_dat_reaction_to_pathway, fscore_pf_reactions, \
        precision_pathways_dat_reaction_to_pathway, precision_pf_reactions, \
        size_pathways_dat_reaction_to_pathway, size_pf_reactions = \
            process_main(pgdb, pre_savi_folder, pathologic_folder, release_folder)

        print(pgdb, len(pre_savi_pathways_dat_reaction_to_pathway), len(reverse_dict(pre_savi_pathways_dat_reaction_to_pathway)), len(pre_savi_pf_reactions),
              len(release_pathways_dat_reaction_to_pathway), len(reverse_dict(release_pathways_dat_reaction_to_pathway)), len(""),
              len(fscore_pathways_dat_reaction_to_pathway), len(reverse_dict(fscore_pathways_dat_reaction_to_pathway)), len(fscore_pf_reactions),
              len(precision_pathways_dat_reaction_to_pathway), len(reverse_dict(precision_pathways_dat_reaction_to_pathway)), len(precision_pf_reactions),
              len(size_pathways_dat_reaction_to_pathway), len(reverse_dict(size_pathways_dat_reaction_to_pathway)), len(size_pf_reactions))