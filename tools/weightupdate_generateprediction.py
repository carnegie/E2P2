import os
from argparse import ArgumentParser


def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))


def generate_prediction_file_for_assessment(fasta_input, e2p2_out, output_folder):
	prediction_path = os.path.join(output_folder, os.path.basename(e2p2_out) + '.prediction.txt')
	training_label_path = os.path.join(output_folder, os.path.basename(e2p2_out) + '.training_label.txt')
	training_instance_path = os.path.join(output_folder, os.path.basename(e2p2_out) + '.training_instance.txt')
	actual_dict = {}
	all_seq_ids = set()
	all_classes = set()
	with open(fasta_input, 'r') as fi:
		for name, seq in read_fasta(fi):
			try:
				info = name[1:].split('|')
				seq_id = info[0].strip()
				all_seq_ids.add(seq_id)
				actual_classes = [c.strip() for c in info[1:] if len(c.strip()) > 0]
				actual_dict.setdefault(seq_id, actual_classes)
				# all_classes = all_classes.union(set(actual_classes))
			except IndexError:
				print('Exclude Line:', name, seq)

	with open(e2p2_out, 'r') as eo, open(prediction_path, 'w') as pp, open(training_label_path, 'w') as tlp, open(training_instance_path, 'w') as tip:
		for line in eo:
			if line.strip().startswith('#'):
				continue
			info = line.split('\t')
			try:
				seq_id = info[0].strip()
				all_seq_ids.add(seq_id)
				predicted_classes = [c.strip() for c in info[1].split('|') if len(c.strip()) > 0]
				all_classes = all_classes.union(set(predicted_classes))
				# for c in predicted_classes:
				# 	tlp.write(c.strip() + '\n')
				try:
					write_line = ""
					actual_classes = actual_dict[seq_id]
					if len(actual_classes) == 0:
						write_line += seq_id + '|NA\t'
					else:
						write_line += seq_id + '|' + '|'.join(sorted(set(actual_classes))) + '\t'
					if len(predicted_classes) == 0:
						write_line += seq_id + '|NA'
					else:
						write_line += seq_id + '|' + '|'.join(sorted(set(predicted_classes)))
					pp.write(write_line + '\n')
				except KeyError:
					print("ID not found", seq_id)
			except IndexError:
				print("Exclude Line:", line)
		for seq in sorted(set(actual_dict.keys()) - set(all_seq_ids)):
			write_line = ""
			actual_classes = actual_dict[seq]
			if len(actual_classes) == 0:
				write_line += seq + '|NA\t'
			else:
				write_line += seq + '|' + '|'.join(sorted(set(actual_classes))) + '\t'
			write_line += seq + '|NA'
			pp.write(write_line + '\n')
		for seqs in sorted(all_seq_ids):
			tip.write(seqs + '\n')
		# for seqs in sorted(actual_dict.keys()):
		# 	tip.write(seqs + '\n')
		all_classes.add('NA')
		for classes in sorted(all_classes):
			tlp.write(classes + '\n')


if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('fasta_input')
	parser.add_argument('classifier_output')
	parser.add_argument('output_folder')

	args = parser.parse_args()
	generate_prediction_file_for_assessment(args.fasta_input, args.classifier_output, args.output_folder)
	# generate_prediction_file_for_assessment(
	# 	'/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-pipeline/UpdateWeights/input/rpsd-4.0-20190108.ef-hold.lst.fa',
	# 	'/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-pipeline/UpdateWeights/classifiers/rpsd-4.0-20190108.ef-hold.lst.fa.e2p2v3.PRIAM',
	# 	'/Users/bxue/Documents/Carnegie/PMNProject/update-rpsd-pipeline/UpdateWeights/predictions')
