import sys, textwrap

out_fasta = sys.argv[1]

in_tsvs = sys.argv[2:]
med_confidence_sets=[]
with open(out_fasta,'w') as out:
	for tsv in in_tsvs:
		for line in open(tsv).readlines()[1:]:
			entry = line.strip().split('\t')
			five_prime=entry[0]
			three_prime=entry[1]
			event=entry[8]
			confidence=entry[16]
			cDNA=entry[20]
			#print(cDNA)
			breakpoint_pos=cDNA.index('|')
			out.write('>fusion.[{}]-[{}].{}.{}.{}\n'.format(five_prime,three_prime,event,confidence,breakpoint_pos))
			cDNA_sanitized=cDNA.replace('|','').replace('_','').replace('.','')
			out.write('{}\n'.format('\n'.join(textwrap.wrap(cDNA_sanitized, 80))))
