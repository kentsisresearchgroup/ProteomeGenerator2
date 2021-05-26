import sys

blast_infile = sys.argv[1]
proteome_infile = sys.argv[2]

all_blast_hits = set()
with open(blast_infile) as blast:
	for line in blast:
		mstrg = line.split('\t')[0]
		all_blast_hits.add(mstrg)

"""
non_canonicals = []
with open(proteome_infile) as proteome:
	for line in proteome:
		if line.startswith('>'):
			#print(line)
			mstrg = line.split(' ')[1]
			if mstrg not in all_blast_hits:
				non_canonicals.append(mstrg)

for mstrg in non_canonicals:
	print(mstrg)
"""

with open(proteome_infile) as proteome:
	for line in proteome:
		if line.startswith('>'):
			non_canonical_sequence = False
			mstrg = line[1:].strip()
			if mstrg not in all_blast_hits:
				print(line.strip())
				non_canonical_sequence = True
		else:
			if non_canonical_sequence:
				print(line.strip())
				

