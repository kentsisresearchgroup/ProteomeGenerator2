import sys

peptides_txt = sys.argv[1]
peptides_fa = sys.argv[2]

peps = open(peptides_txt).readlines()
peps = [x.split('\t') for x in peps[1:]]

for p in peps:
	open(peptides_fa,'a').write('>{}\n{}\n'.format(p[0],p[0]))
