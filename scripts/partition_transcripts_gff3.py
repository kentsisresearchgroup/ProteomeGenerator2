import sys,re
in_gff3 = sys.argv[1]
out_original = sys.argv[2]
out_expanded = sys.argv[3]

with open(in_gff3,'r') as f:
	for line in f:
		if re.search(r'ENST[0-9]+.[0-9]+_[0-9]+',line) is not None:
			open(out_expanded,'a').write(line)
		else:
			open(out_original,'a').write(line)
