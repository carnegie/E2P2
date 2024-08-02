import os

from src.e2p2.classifiers.blast import BLAST


def blast_pred():
    partitions = ['f0', 'f1', 'f2', 'f3', 'f4', 'hold']
    for fold in partitions:
        output = ('/Users/bxue/Documents/Carnegie/PMNProject/RPSDv5.0/E2P2/blastp/' 
                  'rpsd.v5.2.ef.fasta.rpsd.v5.2.ef.ids.txt-' + fold + '.lst.subset.fa.blastp.out')
        b = BLAST('112223', '/Users/bxue/Documents/Carnegie/PMNProject/RPSDv4.2/weights/blast')
        b.read_classifier_result(output)

        file_name, file_ext = os.path.splitext(os.path.basename(output))

        predict = os.path.join('/Users/bxue/Documents/Carnegie/PMNProject/RPSDv5.0/E2P2/predict',
                               'rpsd.v5.2.ef.ids.txt-' + fold + '.blastp.e2p2')
        with open(predict, 'w') as op:
            for query in b.res:
                cls_of_query = set()
                for cls in b.res[query]:
                    cls_of_query.add(cls.name)
                op.write(query + '\t' + '|'.join(sorted(cls_of_query)) + '\n')

