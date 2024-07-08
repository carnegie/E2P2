import os

from argparse import ArgumentParser


def measure_output_helper(fp):
	header, entries = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith("#Performance of each"):
			if header: yield (header, entries)
			header, entries = line, []
		elif not line.startswith('#'):
			entries.append(line)
	if header: yield (header, entries)


def read_ef_f_score(file_path, ef_classes_dict):
	with open(file_path, 'r') as fp:
		for header, entries in measure_output_helper(fp):
			if "label" in header:
				for entry in [e for e in entries if len(e) > 0]:
					info = [e.strip() for e in entry.split('\t')]
					# print(info)
					try:
						ef_classes_dict[info[0]].append(float(info[-1]))
					except KeyError:
						ef_classes_dict.setdefault(info[0], [float(info[-1])])


def calculate_weight(ef_classes_dict, output_path):
	with open(output_path, 'w') as op:
		for ef_class in sorted(ef_classes_dict.keys()):
			weight = sum(ef_classes_dict[ef_class])/len(ef_classes_dict[ef_class])
			# weight = sum(ef_classes_dict[ef_class]) / float(6)
			print(ef_class, ef_classes_dict[ef_class], weight)
			op.write(ef_class + "\t" + "{0:.3f}".format(weight) + "\n")


if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('input_folder')
	parser.add_argument('output_path')

	args = parser.parse_args()
	ef_classes_dict = {}
	for f in os.listdir(args.input_folder):
		file_path = os.path.join(args.input_folder, f)
		if os.path.isfile(file_path) and not f.startswith('.') and f.endswith('.output.txt'):
			read_ef_f_score(file_path, ef_classes_dict)
	calculate_weight(ef_classes_dict, args.output_path)
