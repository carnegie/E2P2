import os
from argparse import ArgumentParser


def reade2p2(e2p2_output):
	e2p2_id, e2p2_attr = "", ""
	for line in e2p2_output:
		if line.startswith("ID\t"):
			if len(e2p2_id) > 0:
				yield e2p2_id, e2p2_attr
			e2p2_id, e2p2_attr = line, ""
		else:
			e2p2_attr += line
	if len(e2p2_id) > 0:
		yield e2p2_id, e2p2_attr


def getidsfrome2p2(file_path):
	pf_dict = {}
	with open(file_path, 'r') as fp:
		for e2p2_id, e2p2_attr in reade2p2(fp):
			if e2p2_id is not None and len(e2p2_id) > 0:
				attrs_dict = pf_dict.setdefault(e2p2_id.split('\t')[1].strip(), {})
				# id_set.add(e2p2_id.rstrip())
				attrs_info = [i.split('\t') for i in e2p2_attr.split('\n') if len(i.strip()) > 0 and i.strip() != '//']
				for attr in attrs_info:
					try:
						attrs_dict[attr[0].strip()].append(attr[1].strip())
					except KeyError:
						attrs_dict.setdefault(attr[0].strip(), [attr[1].strip()])
	return pf_dict


def compareids(file_path_1, file_path_2, out_p1, out_p2):
	pf_dict_1 = getidsfrome2p2(file_path_1)
	file_name_1 = os.path.basename(file_path_1)
	pf_dict_2 = getidsfrome2p2(file_path_2)
	file_name_2 = os.path.basename(file_path_2)

	out_p1.write('Seq\t' + file_name_1 + '\t' + file_name_2 + '\n')
	seq_in1_not2 = sorted(list(set(pf_dict_1.keys()) - set(pf_dict_2.keys())))
	seq_in2_not1 = sorted(list(set(pf_dict_2.keys()) - set(pf_dict_1.keys())))
	seq_in_both = sorted(list(set(pf_dict_1.keys()) & set(pf_dict_2.keys())))

	for seq in sorted(set(pf_dict_1.keys()).union(set(pf_dict_2.keys()))):
		if seq in seq_in1_not2:
			out_p1.write(seq + '\t' + 'Y' + '\t' + 'N' + '\n')
		elif seq in seq_in2_not1:
			out_p1.write(seq + '\t' + 'N' + '\t' + 'Y' + '\n')
		elif seq in seq_in_both:
			out_p1.write(seq + '\t' + 'Y' + '\t' + 'Y' + '\n')
		else:
			out_p1.write(seq + '\t' + 'N' + '\t' + 'N' + '\n')

	# if len(seq_in1_not2) > 0:
	# 	print('Seq ID in file 1 not file 2', seq_in1_not2)
	# if len(seq_in2_not1) > 0:
	# 	print('Seq ID in file 2 not file 1', seq_in2_not1)

	out_p2.write('Seq\tAttr\tOnly In ' + file_name_1 + '\tOnly In ' + file_name_2 + '\tIn Both\n')
	for seq_id in seq_in_both:
		attrs_dict_1 = pf_dict_1[seq_id]
		attrs_dict_2 = pf_dict_2[seq_id]

		# attrs_in1_not2 = sorted(list(set(attrs_dict_1.keys()) - set(attrs_dict_2.keys())))
		# attrs_in2_not1 = sorted(list(set(attrs_dict_2.keys()) - set(attrs_dict_1.keys())))
		attrs_in_both = sorted(list(set(attrs_dict_1.keys()) & set(attrs_dict_2.keys())))

		# if len(attrs_in1_not2) > 0:
		# 	for attr1 in attrs_in1_not2:
		# 		print('For Seq', seq_id, 'Attribute in file 1 not file 2', attrs_in1_not2, attrs_dict_1[attr1])
		# if len(attrs_in2_not1) > 0:
		# 	for attr2 in attrs_in2_not1:
		# 		print('For Seq', seq_id, 'Attribute in file 2 not file 1', attrs_in2_not1, attrs_dict_2[attr2])

		for attr in attrs_in_both:
			if attr != 'NAME' and attr != 'PRODUCT-TYPE':
				value_list_1 = attrs_dict_1[attr]
				value_list_2 = attrs_dict_2[attr]
				vals_in1_not2 = sorted(list(set(value_list_1) - set(value_list_2)))
				vals_in2_not1 = sorted(list(set(value_list_2) - set(value_list_1)))
				vals_in_both = sorted(list(set(value_list_1) & set(value_list_2)))
				# print(seq_id, vals_in_both)

				# if len(vals_in1_not2) > 0:
				# 	print('For Seq', seq_id, 'Attribute', attr, 'Value in file 1 not file 2', vals_in1_not2)
				# if len(vals_in2_not1) > 0:
				# 	print('For Seq', seq_id, 'Attribute', attr, 'Value in file 2 not file 1', vals_in2_not1)

				out_p2.write(seq_id + '\t' + attr + '\t' + ','.join(vals_in1_not2) + '\t' + ','.join(vals_in2_not1) + '\t' + ','.join(vals_in_both) + '\n')


if __name__ == '__main__':
	parser = ArgumentParser()
	parser.add_argument('e2p2_1')
	parser.add_argument('e2p2_2')

	args = parser.parse_args()

	print('File 1', args.e2p2_1)
	print('File 2', args.e2p2_2)
	output_1 = args.e2p2_1 + '.' + os.path.basename(args.e2p2_2) + '.diff_seq.txt'
	output_2 = args.e2p2_1 + '.' + os.path.basename(args.e2p2_2) + '.diff_val.txt'
	with open(output_1, 'w') as op1, open(output_2, 'w') as op2:
		compareids(args.e2p2_1, args.e2p2_2, op1, op2)
