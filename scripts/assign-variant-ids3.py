import sys
from variantclassifier import determine_class_from_line

prev_chrom = ''
prev_pos = ''

counter = 0

for line in sys.stdin:
	if line.startswith('##'):
		print(line[:-1])
		continue
	if line.startswith('#'):
		print('##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs per ALT allele.\">')
		print(line[:-1])
		continue
	fields = line.split()
	chrom = fields[0]
	start = fields[1]
	vartype = determine_class_from_line(line)
	alleles = [fields[3]] + [a for a in fields[4].split(',')]
	length = max(1, max([len(a) for a in alleles])-1)
	if (chrom == prev_chrom):
		counter += 1
	else:
		counter = 0
	var_id = str(counter)
	fields[2] = var_id
	prev_chrom = chrom
	prev_pos = start
	prev_len = length
	fields[7] = 'ID=' + var_id
	print('\t'.join(fields))
