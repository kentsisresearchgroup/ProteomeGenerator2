import sys
peptides_txt = sys.argv[1]
#blast_outfmt6 = sys.argv[2]
ref_db = sys.argv[2]
novel_peps_out = sys.argv[3]
novel_mstrg_map_out = sys.argv[4]

peps = open(peptides_txt).readlines()
peps = [x.split('\t') for x in peps[1:]]
peps_set = set()
for p in peps:
	peps_set.add(p[0])
"""
blastHits = open(blast_outfmt6).readlines()
blastHits=[x.split('\t') for x in blastHits]
blast_set = set()
for p in blastHits:
	blast_set.add(p[0])
"""
seqs=[]
seq=''
#with open('/data/kentsis/indexes/UniProt/UP_sprot_human_canonical+isoforms+cRAP.fasta') as uniprot:
with open(ref_db) as ref:
	for line in ref:
		line = line.strip()
		if not line:
			continue
		if line.startswith('>'):
			if seq is not None:
				seqs.append(seq)
				seq=''
			continue
		seq=seq+line

nonhits = [pep for pep in peps_set if sum([1 for seq in seqs if pep in seq])==0]
nonhits_full = [x for x in peps if x[0] in nonhits]

rev_index = -1
for p in peps:
	for element in p:
		if 'REV_' in element:
			rev_index = p.index(element)
	if rev_index != -1: break

nonhits_noREV_noCON = [x for x in nonhits_full if 'REV_' not in x[rev_index] and 'CON_' not in x[rev_index-1]]

for p in nonhits_noREV_noCON:
	open(novel_peps_out,'a').write('{}\n'.format(p[0]))
	open(novel_mstrg_map_out,'a').write('{}\t{}\n'.format(p[0],p[rev_index-1]))


