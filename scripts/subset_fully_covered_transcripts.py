import sys,re

output_filename = sys.argv[1]
annotation_gtf = open(sys.argv[2]).readlines()

covered_transcripts_set = set()
first_pass = True
for covered_transcripts_gtf in sys.argv[3:]:
	covered_transcripts = open(covered_transcripts_gtf).readlines()
	updated_covered_transcripts_set = set()
	prev_line = ''
	
	for line in covered_transcripts:
		if 'ID=' in line:
			#print(line)
			if prev_line and line.startswith(';'): line = prev_line + line
			#print(line)
			transcript=re.search(r'EN[A-Z]*ST[0-9]+.[0-9]+',line).group()
			if first_pass:
				#print(transcript)
				updated_covered_transcripts_set.add(transcript)
			else:
				if transcript in covered_transcripts_set:
					updated_covered_transcripts_set.add(transcript)
		prev_line = line
	covered_transcripts_set = updated_covered_transcripts_set
	#print(len(covered_transcripts_set))
	first_pass = False
for line in annotation_gtf:
	match = re.search(r'EN[A-Z]*ST[0-9]+.[0-9]+',line)
	if match:
		enst = match.group()
		if enst in covered_transcripts_set:
			open(output_filename,'a').write(line)
