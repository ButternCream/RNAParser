import re, os, argparse, shutil
# Chart reference for amino acids
amino_acid_chart = {
	"AAA": "K",
	"AAC": "N",
	"AAG": "K",
	"AAU": "N",
	"ACA": "T",
	"ACC": "T",
	"ACG": "T",
	"ACU": "T",
	"AGA": "R",
	"AGC": "S",
	"AGG": "R",
	"AGU": "S",
	"AUA": "I",
	"AUC": "I",
	"AUG": "M",
	"AUU": "I",
	"CAA": "Q",
	"CAC": "H",
	"CAG": "Q",
	"CAU": "H",
	"CCA": "P",
	"CCC": "P",
	"CCG": "P",
	"CCU": "P",
	"CGA": "R",
	"CGC": "R",
	"CGG": "R",
	"CGU": "R",
	"CUA": "L",
	"CUC": "L",
	"CUG": "L",
	"CUU": "L",
	"GAA": "E",
	"GAC": "D",
	"GAG": "E",
	"GAU": "D",
	"GCA": "A",
	"GCC": "A",
	"GCG": "A",
	"GCU": "A",
	"GGA": "G",
	"GGC": "G",
	"GGG": "G",
	"GGU": "G",
	"GUA": "V",
	"GUC": "V",
	"GUG": "V",
	"GUU": "V",
	"UAA": "STOP",
	"UAC": "Y",
	"UAG": "STOP",
	"UAU": "Y",
	"UCA": "S",
	"UCC": "S",
	"UCG": "S",
	"UCU": "S",
	"UGA": "STOP",
	"UGC": "C",
	"UGG": "W",
	"UGU": "C",
	"UUA": "L",
	"UUC": "F",
	"UUG": "L",
	"UUU": "F"
}

# Create arguments
parser = argparse.ArgumentParser(description='Converts files which contain an RNA sequence and outputs the corresponding peptides to a file')
parser.add_argument('-e', '--extension', dest='extension',default='pep', nargs=1, help='Specify the file extension for the outputted files')
parser.add_argument('-o', '--output', dest='output', type=str, default='Peptides', nargs=1, help='Specify the output directory')
parser.add_argument('-f', '--files', dest='files', type=argparse.FileType('r'),nargs='+', help='Specify the files to read')
parser.add_argument('-c', '--clean', dest='clean', nargs='+',help='Clean the specified directories')

class Peptide:
	def __init__(self, acid_chain, seq):
		"""
		Constructor
		:param acid_chain: A string chain of acids, e.g MKNTT (from acid chart table)
		:param seq: The corresponding sequence in the RNA string
		"""
		self.acid_chain = acid_chain
		self.seq = seq

def translate(seq, start_index=0):
	"""
	Translate a subset of the sequence provided into list of amino acids
	:param seq: The RNA sequence to parse as a string
	:param start_index: Index to start at in the sequence
	:return: list of amino acids found in parsed sequence, where it started, where it ended, if it failed (i.e never found 'STOP' amino acid)
	TODO: Simply pass in the substring instead of a new starting index
	"""
	b = start_index
	found_start = False
	found_finish = False
	end_index = 0
	aa_list = list()
	while b < len(seq) and not found_finish:
		if not found_start and len(seq[b:b+3]) == 3 and amino_acid_chart[seq[b:b+3]] == "M":
			start_index = b
			found_start = True
		if len(seq[b:b+3]) == 3 and amino_acid_chart[seq[b:b+3]] == "STOP" and found_start:
			end_index = b
			found_finish = True
		if found_start and len(seq[b:b + 3]) == 3 and not found_finish:
			aa_list.append(amino_acid_chart[seq[b:b + 3]])
		b += 3

	if not found_finish:
		return list(), start_index, b , True
	return aa_list, start_index, end_index, False


def get_peptides(seq):
	"""
	Get list of all peptides found in the RNA sequence
	:param seq: The RNA sequence as a string
	:return: List of peptides, ending index of the last successful translation (AUGUUUUGAAUGUUUUUU would return 8 as the last successful index)
	"""
	last = 0
	peptide_list = list()
	l,s,e,f = translate(seq)
	if not f:
		peptide_list.append(Peptide(''.join(l), seq[s:e+3]))
	while s < len(seq) and not f:
		last = e
		'''
		l = list
		s = start index
		e = end index
		f = fail
		'''
		l, s, e, f = translate(seq, s+3)
		if not f:
			peptide_list.append(Peptide(''.join(l), seq[s:e+3]))
		s = e
	return peptide_list, last

def print_peptides(p_list):
	"""
	Print out the peptides in the given list
	:param p_list: List of peptides
	:return:
	"""
	for p in p_list:
		print("Peptide:", p.acid_chain)
		print("Sequence:", p.seq)
		print()
	print("Total Peptides:", len(p_list))

def convert(rRNA):
	"""
	Converts RNA sequence with T's to have U's, if it doesn't already
	:param rRNA: The RNA sequence to convert as a string
	:return: Converted RNA sequence or -1 if it contains invalid characters
	"""
	if re.match("[AGCU]", rRNA):
		return rRNA
	elif re.match("[^AGCT]", rRNA):
		return -1
	new = rRNA.replace('T', 'U')
	new = re.sub(r"\s+", '', new)
	return new

def create_dir(d):
	if not os.path.exists(d):
		os.makedirs(d)

def clean(d_list):
	[shutil.rmtree(d) for d in d_list if os.path.exists(d)]

def main():
	args = parser.parse_args()
	if args.clean is not None:
		clean(args.clean)
		return
	if args.files is not None:
		if type(args.output) is list: output_dir = args.output[0]
		else: output_dir = args.output
		create_dir(output_dir)
		for f in args.files:
			data=f.read().replace('\n', '')
			f.close()
			seq = convert(data)
			if seq != -1:
				p_list, last = get_peptides(seq)
				new_file = str(os.path.splitext(f.name)[0].split('/')[-1]) + '.{}'.format(args.extension)
				new_file = os.path.join(output_dir, new_file)
				with open(new_file, 'w') as new_f:
					new_f.write('Original File: {0}\n\n'.format(f.name))
					for p in p_list:
						text = 'Peptide: {0}\nSubset of RNA Sequence: {1}\n'.format(p.acid_chain, p.seq)
						new_f.write(text)
					new_f.write("\nParsed through %.2f%% of RNA Sequence\n" % (last / len(seq) * 100))
					new_f.close()
			else:
				print('Invalid character found in {0}'.format(f.name))

main()