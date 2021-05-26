import sys, re
from Bio import pairwise2

novel_proteins_fasta = sys.argv[1]

"""
nonmutational_novel_peptides_file = sys.argv[1]
peptide_mstrgs_infile = sys.argv[2]
"""

"""
proteome_fa_infile = sys.argv[2]
"""

ref_fa_infile = sys.argv[2]
mstrg_ref_blast_infile = sys.argv[3]

genome_annotation_infiles = sys.argv[4:]

"""
nonmutational_novel_peps = [x.strip() for x in open(nonmutational_novel_peptides_file).readlines()]
"""

"""
peptide_mstrgs_dict = dict()
for line in open(peptide_mstrgs_infile).readlines():
	line = line.strip()
	peptide = line.split('\t')[0]
	mstrgs = line.split('\t')[1]

	peptide_mstrgs_dict[peptide] = mstrgs.split(';')
"""

novel_proteome_dict = dict()
novel_proteome_mstrgs = []
seq=''
header=''
with open(novel_proteins_fasta) as f:
	for line in f:
		line = line.strip()
		if not line: continue
		if line.startswith('>'):
			if seq is not None:
				novel_proteome_dict[header]=seq
			seq=''
			header=line[4:-1].split(' ')[0]
			novel_proteome_mstrgs.append(header)
		else:
			seq=seq+line
novel_proteome_dict[header]=seq


"""
proteome_dict=dict()
seq=''
header=''
with open(proteome_fa_infile) as f:
        for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith('>'):
                        if seq is not None:
                                proteome_dict[header]=seq
                        seq=''
                        header=line[4:-1].split(' ')[0]
                else:
                        seq=seq+line
proteome_dict[header]=seq
"""


ref_dict=dict()
seq=''
header=''
with open(ref_fa_infile) as ref:
        for line in ref:
                line = line.strip()
                if not line: continue
                if line.startswith('>'):
                        if seq is not None:
                                ref_dict[header]=seq
                                seq=''
                        header=re.search(r'ENST[0-9]+',line).group()
                else:
                        seq=seq+line
ref_dict[header]=seq

header=''
mstrg_ref_blast_dict = dict()
for line in open(mstrg_ref_blast_infile).readlines():
	entry = line.strip().split('\t')
	mstrg = entry[0][3:-1]
	alignment_len = int(entry[3])
	mismatches = int(entry[4])
	if mismatches > 0.10*alignment_len: continue
	enst = re.search(r'ENST[0-9]+',entry[1]).group()

	if mstrg in mstrg_ref_blast_dict:
		lst = mstrg_ref_blast_dict[mstrg]
		lst.append(enst)
		mstrg_ref_blast_dict[mstrg] = lst
	else:
		mstrg_ref_blast_dict[mstrg] = [enst]

"""
mstrg_genomic_dict = dict()
with open(peptide_mstrgs_infile) as f:
	for line in f:
		if line.split('\t')[0] in nonmutational_novel_peps:
			mstrgs = line.strip().split('\t')[1].split(';')
			for mstrg in mstrgs:
				if 'chr' not in mstrg or 'mm10' in mstrg: continue
				#print(mstrg)
				chrom = re.search(r'[12]_chr[0-9XY]+',mstrg).group()

				coords = re.search(r'[0-9]+-[0-9]+',mstrg).group()
				l_coord = int(coords.split('-')[0])
				r_coord = int(coords.split('-')[1])

				strand = re.search(r'\([+-]\)',mstrg).group()[1]

				if chrom not in mstrg_genomic_dict:
					mstrg_genomic_dict[chrom] = [(l_coord, r_coord, strand, mstrg)]
				else:
					coords_list = mstrg_genomic_dict[chrom]
					coords_list.append((l_coord, r_coord, strand, mstrg))
					mstrg_genomic_dict[chrom] = coords_list

	for chrom in mstrg_genomic_dict:
		coords_list = mstrg_genomic_dict[chrom]
		coords_list.sort(key=lambda tup:tup[0]) # sort by l_coord
		mstrg_genomic_dict[chrom] = coords_list
"""

mstrg_genomic_dict = dict()
for mstrg in novel_proteome_mstrgs:
	if 'chr' not in mstrg or 'mm10' in mstrg: continue
	#print(mstrg)

	match = re.search(r'[12]_chr[0-9XYM]+',mstrg)
	if not match: continue
	chrom = match.group()
	#chrom = re.search(r'[12]_chr[0-9XYM]+',mstrg).group()

	coords = re.search(r'[0-9]+-[0-9]+',mstrg).group()
	l_coord = int(coords.split('-')[0])
	r_coord = int(coords.split('-')[1])

	strand = re.search(r'\([+-]\)',mstrg).group()[1]

	if chrom not in mstrg_genomic_dict:
		mstrg_genomic_dict[chrom] = [(l_coord, r_coord, strand, mstrg)]
	else:
		coords_list = mstrg_genomic_dict[chrom]
		coords_list.append((l_coord, r_coord, strand, mstrg))
		mstrg_genomic_dict[chrom] = coords_list

	for chrom in mstrg_genomic_dict:
		coords_list = mstrg_genomic_dict[chrom]
		coords_list.sort(key=lambda tup:tup[0]) # sort by l_coord
		mstrg_genomic_dict[chrom] = coords_list

lncRNA_set=set()
pseudogene_set=set()
anti_sense_set=set()
upstreamARF_set=set()
fivePrimeUtrORF_set=set()
threePrimeUtrORF_set=set()
denovo_transcript_set=set()
for infile in genome_annotation_infiles:
	with open(infile) as gtf:
		gene_of_interest = False
		for line in gtf:

			m = re.search(r'H([12]).gtf',infile)
			haplo = '' if not m else m.group(1)
			entry = line.split('\t')
			feature_type = entry[2]
	
			if (not gene_of_interest) and feature_type != 'gene': continue	
		
			chrom = entry[0]
			l_coord_gtf = int(entry[3])
			r_coord_gtf = int(entry[4])
			strand_gtf = entry[6]

			info = entry[8]
			
			if feature_type == 'gene': 

				# resolve features of previous gene
				if gene_of_interest:
					coords_list = mstrg_genomic_dict[current_gene['chrom']]
					for l_coord_mstrg, r_coord_mstrg, strand_mstrg, mstrg in coords_list:
						if haplo and '{}_chr'.format(haplo) not in mstrg: continue
						if current_gene['l_bound'] <= l_coord_mstrg and current_gene['r_bound'] >= r_coord_mstrg:
							for t in current_gene['transcripts']:
								if strand_mstrg != current_gene['strand']: 
									anti_sense_set.add(mstrg)
									continue
								if t not in current_gene['CDS']: continue # temporary workaround
								CDS_exons = current_gene['CDS'][t]
								if 'UTR' not in current_gene or t not in current_gene['UTR']: continue
								UTR_exons = current_gene['UTR'][t]
								if strand_mstrg=='+':
									CDS_exons.sort(key=lambda tup:tup[0])
									UTR_exons.sort(key=lambda tup:tup[0])

									#print(UTR_exons)

									fiveP=[x for x in UTR_exons if x[1]<=CDS_exons[0][0]]
									threeP=[x for x in UTR_exons if x[0]>=CDS_exons[-1][1]]
									#print(fiveP)
									#print(threeP)
									#print(t)
									#print('left_bound = {}'.format(l_coord_mstrg))
									#print('right_bound = {}'.format(r_coord_mstrg))
									#print('3prime left bound = {}'.format(threeP[0][0]))
									#print('3prime right bound = {}'.format(threeP[-1][1]))
									if fiveP and fiveP[0][0] < l_coord_mstrg and CDS_exons[0][0] < r_coord_mstrg and CDS_exons[0][1] > r_coord_mstrg:
										upstreamARF_set.add(mstrg)
									elif fiveP and fiveP[0][0] < l_coord_mstrg and fiveP[-1][1] > r_coord_mstrg:
										fivePrimeUtrORF_set.add(mstrg)
									elif threeP and threeP[0][0] < l_coord_mstrg and threeP[-1][1] > r_coord_mstrg:
										threePrimeUtrORF_set.add(mstrg)
									else:
										denovo_transcript_set.add(mstrg)
							
								elif strand_mstrg=='-':
									CDS_exons.sort(key=lambda tup:tup[1],reverse=True)
									UTR_exons.sort(key=lambda tup:tup[1],reverse=True)

									fiveP=[x for x in UTR_exons if x[0]>=CDS_exons[0][1]]
									threeP=[x for x in UTR_exons if x[1]<=CDS_exons[-1][0]]

									if fiveP and fiveP[0][1] > r_coord_mstrg and CDS_exons[0][1] > l_coord_mstrg and CDS_exons[0][0] < l_coord_mstrg:
										upstreamARF_set.add(mstrg)
									elif fiveP and fiveP[0][1] > r_coord_mstrg and fiveP[-1][0] < l_coord_mstrg:
										fivePrimeUtrORF_set.add(mstrg)
									elif threeP and threeP[0][1] > r_coord_mstrg and threeP[-1][0] < l_coord_mstrg:
										threePrimeUtrORF_set.add(mstrg)
									else:
										denovo_transcript_set.add(mstrg)
				
				gene_of_interest = False
				current_gene=dict()
				if chrom not in mstrg_genomic_dict: continue
				coords_list = mstrg_genomic_dict[chrom]
				for l_coord_mstrg, r_coord_mstrg, strand_mstrg, mstrg in coords_list:
					if l_coord_gtf <= l_coord_mstrg and r_coord_gtf >= r_coord_mstrg:
						if 'lncRNA' in info: lncRNA_set.add(mstrg)
						elif 'pseudogene' in info: pseudogene_set.add(mstrg)	
						else:
							gene_of_interest = True
							current_gene['chrom'] = chrom
							current_gene['l_bound'] = l_coord_gtf
							current_gene['r_bound'] = r_coord_gtf
							current_gene['strand'] = strand_gtf
							current_gene['transcripts'] = set()

			elif feature_type == 'UTR' or feature_type == 'CDS':
				enst = re.search(r'ENST[0-9]+',info).group()
				s = current_gene['transcripts']
				s.add(enst)
				current_gene['transcripts'] = s

				if feature_type not in current_gene: current_gene[feature_type] = {enst: [(l_coord_gtf,r_coord_gtf)]}
				else: 
					feature_dict=current_gene[feature_type] 
					if enst in feature_dict:
						lst = feature_dict[enst]
						lst.append((l_coord_gtf,r_coord_gtf))
						feature_dict[enst]=lst
					else:
						feature_dict[enst] = [(l_coord_gtf,r_coord_gtf)]
					current_gene[feature_type] = feature_dict


print('upstreamARFs: {}'.format(upstreamARF_set))
print('5prime orfs: {}'.format(fivePrimeUtrORF_set))
print('3prime orfs: {}'.format(threePrimeUtrORF_set))
print('antisense orfs: {}'.format(anti_sense_set))
print('denovo transcripts: {}'.format(denovo_transcript_set))
print('lncRNAs: {}'.format(lncRNA_set))
print('psudogenes: {}'.format(pseudogene_set))

"""
for line in open(nonmutational_novel_peptides_file).readlines():
	pep = line.strip()

	variant_type=set()
	top_alignment_formatted = '' 
	for mstrg in peptide_mstrgs_dict[pep]:

		if mstrg in mstrg_ref_blast_dict:

			for enst in mstrg_ref_blast_dict[mstrg]:
				mstrg_seq = proteome_dict[mstrg]
				ref_seq = ref_dict[enst]

				alignments=pairwise2.align.localmd(pep,ref_seq,4,0.1,-3,-2.5,0,0, penalize_end_gaps=(False,True))
				top_alignment_formatted = pairwise2.format_alignment(*alignments[0])
				top_align_lst = [x for x in top_alignment_formatted.split('\n')]
				mstrg_align = top_align_lst[0].split(' ')[-1]
				ref_align = top_align_lst[2].split(' ')[-1]
				comps = top_align_lst[1][-1*len(mstrg_align):]
				
				if len(pep) == len(comps) and sum([1 for x in comps if x=='|']) == len(comps): # perfect match
					variant_type.add('canonical_peptide')
					#break

				elif len(pep) == len(comps) and sum([1 for x in comps if x=='.'])==1 and sum([1 for x in comps if x=='|'])==len(comps)-1:
					variant_type.add('offset_exon_boundary or missense')
					#break
				elif int(0.4*len(comps))*'|' in comps and (len(comps)>=len(pep) and sum([1 for x in comps if x=='|'])+sum([1 for x in comps if x==' ']) == len(comps)):
					variant_type.add('offset_exon_boundary or indel')
					#variant_type.add('likely_indel')
					#break
				elif mstrg in pseudogene_set: 
					variant_type.add('pseudogene')
					#break

				else: # nonsensical alignment
					a=pairwise2.align.localmd(mstrg_seq,ref_seq,4,0.1,-3,-2.5,0,0, penalize_end_gaps=(False,True))
					
					a_formatted=pairwise2.format_alignment(*a[0])
					#print(a_formatted)
					#print(a)
					a_lst = [x for x in a_formatted.split('\n')]
					m_align = a_lst[0].split(' ')[-1]
					r_align = a_lst[2].split(' ')[-1]
					
					#m_start_pos = int(a_lst[0].split(' ')[-2])
					#r_start_pos = int(a_lst[2].split(' ')[-2])
					cmps = a_lst[1][-1*len(r_align):]
					#print(a)
					#print(ref_seq)
					if mstrg_seq.index(pep) < a[0][3] and (a[0][4]-a[0][3])+1 >= 0.8*len(ref_seq):
					#if r_start_pos == 1 and mstrg_seq.index(pep) < m_start_pos:
						#print(a_formatted)
						variant_type.add('upstream_start')
						#break
						

					
					
		if True:
		#else: # blast MISS
		
			# use 
			if mstrg in lncRNA_set: variant_type.add('lncRNA')
			elif mstrg in pseudogene_set: variant_type.add('pseudogene')
			elif mstrg in upstreamARF_set: variant_type.add('5prime_alternate_reading_frame')
			elif mstrg in fivePrimeUtrORF_set: variant_type.add('5prime_UTR_ORF')
			elif mstrg in threePrimeUtrORF_set: variant_type.add('3prime_UTR_ORF')
			elif mstrg in anti_sense_set: variant_type.add('anti-sense_ORF') 
			elif mstrg in denovo_transcript_set: variant_type.add('denovo_transcript')

		#if variant_type != '': break
	if len(variant_type) == 0: variant_type.add('noncanonical_transcript')
	#if top_alignment_formatted: print(top_alignment_formatted)
	print("{}\t{}\t{}".format(pep,','.join(list(variant_type)),peptide_mstrgs_dict[pep]))
"""
