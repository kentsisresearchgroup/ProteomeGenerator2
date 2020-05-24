import sys

fasta= open(sys.argv[1]).readlines()

for line in fasta:
	fixed=line.strip().replace(' ','')
	print(fixed)
