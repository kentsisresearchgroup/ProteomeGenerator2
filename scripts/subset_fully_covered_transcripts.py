import sys,re

output_filename = sys.argv[1]
annotation_gtf = open(sys.argv[2]).readlines()

covered_transcripts_set = set()

for covered_transcripts_gtf in sys.argv[3:]:
	covered_transcripts = open(covered_transcripts_gtf).readlines()
	prev_line = ''
	for line in covered_transcripts:
		if 'ID=' in line:
			#print(line)
			if prev_line and line.startswith(';'): line = prev_line + line
			transcript=re.search(r'ENST[0-9]+.[0-9]+',line).group()
			covered_transcripts_set.add(transcript)
		prev_line = line

for line in annotation_gtf:
	match = re.search(r'ENST[0-9]+.[0-9]+',line)
	if match:
		enst = match.group()
		if enst in covered_transcripts_set:
			open(output_filename,'a').write(line)
