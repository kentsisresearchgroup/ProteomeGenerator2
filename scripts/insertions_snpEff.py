import sys, re, difflib
from Bio import pairwise2

proteome_infile = sys.argv[1]
annovar_infile = sys.argv[2]
blast_outfmt6_infile = sys.argv[3]
ref_fasta_file = sys.argv[4]

novelPep_mstrgs_infile = sys.argv[5]

mutation_mstrg_map_outfile = sys.argv[6]
mutation_MQevidence_map_outfile = sys.argv[7]
MQnovelPep_mutation_map_outfile = sys.argv[8]

blast_dict=dict()
enst_uniprot_dict=dict()
ref_blastHits=set()
with open(blast_outfmt6_infile) as blast:
        for hit in blast:
                mstrg = hit.split('\t')[0][3:-1]
                ref_header = hit.split('\t')[1]
                enst = re.search(r'ENST[0-9]+',ref_header).group()
                enst_uniprot_dict[enst] = ref_header.split('|')[5]
                if enst in blast_dict:
                        lst = blast_dict[enst]
                        lst.append(mstrg)
                        blast_dict[enst] = lst
                else:
                        blast_dict[enst] = [mstrg]
                ref_blastHits.add(enst)

ref_dict=dict()
seq=''
header=''
with open(ref_fasta_file) as ref:
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
seq=''
header=''
proteome_dict=dict()
with open(proteome_infile) as f:
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

mstrg_MQnovelPeps_dict = dict()
with open(novelPep_mstrgs_infile) as f:
        for line in f:
                line = line.strip()
                pep = line.split('\t')[0]
                mstrgs = line.split('\t')[1]

                for m in mstrgs.split(';'):
                        if m not in mstrg_MQnovelPeps_dict:
                                mstrg_MQnovelPeps_dict[m] = [pep]
                        else:
                                pep_lst = mstrg_MQnovelPeps_dict[m]
                                pep_lst.append(pep)
                                mstrg_MQnovelPeps_dict[m] = pep_lst


AA_lookup = dict()
AA_lookup['A'] = 'Ala'
AA_lookup['C'] = 'Cys'
AA_lookup['D'] = 'Asp'
AA_lookup['E'] = 'Glu'
AA_lookup['F'] = 'Phe'
AA_lookup['G'] = 'Gly'
AA_lookup['H'] = 'His'
AA_lookup['I'] = 'Ile'
AA_lookup['K'] = 'Lys'
AA_lookup['L'] = 'Leu'
AA_lookup['M'] = 'Met'
AA_lookup['N'] = 'Asn'
AA_lookup['P'] = 'Pro'
AA_lookup['Q'] = 'Gln'
AA_lookup['R'] = 'Arg'
AA_lookup['S'] = 'Ser'
AA_lookup['T'] = 'Thr'
AA_lookup['V'] = 'Val'
AA_lookup['W'] = 'Trp'
AA_lookup['Y'] = 'Tyr'
AA_lookup['X'] = 'XXX'
AA_lookup['*'] = 'Ter'


revAA_lookup = {v: k for k, v in AA_lookup.items()}


total_insertion = 0
variants_with_expressed_transcripts = 0
successes=0
true_fails=0
true_fail_transcripts=[]
true_fail_set=set()
fails_1_isoform=0
false_fail_transcripts=[]
unsearchables = 0
unsearchable_transcripts=[]
transcripts_missing_AA=set()
non_ref_enst = []
with open(annovar_infile,'r') as f:
    for line in f:
        chrom = line.split('\t')[0]
        successfully_recovered_variant = False
        at_least_1_tr_w_2_isoforms=False
        at_least_1_transcript_without_repeats=False
        variant_has_matched_transcript = False
        match = re.search(r'inframe_insertion',line)
        if match is None: continue
        total_insertion = total_insertion+1
        transcripts = re.findall(r'ENST[0-9]+',line)
        for t in transcripts:
            coding = re.search(r'protein_coding',line.split(t)[1].split(',')[0])
            if not coding: continue    
            if t not in ref_dict:
                print('{}: transcript not in ref protein fasta') 
                continue
            print(t)
            AA_change_full = re.search(r'p\.[A-Z][a-z]+[0-9]+_[A-Z][a-z]+[0-9]+ins([A-Z][a-z]+)+',line.split(t)[1].split(',')[0])
            if AA_change_full is None:
                AA_change_full = re.search(r'p\.([A-Z][a-z][a-z][0-9]+_)?[A-Z][a-z][a-z][0-9]+dup',line.split(t)[1].split(',')[0])
                if AA_change_full is None:
                    print("couldn't find AA change: {}".format(t))
                    transcripts_missing_AA.add(t)
                    bailedOut = True
                    continue
            AA_change_full = AA_change_full.group()
            print(AA_change_full)

            left=''
            right=''
            l_AA=''
            r_AA=''
            insertion=''
            if 'ins' in AA_change_full:
                left = AA_change_full.split('_')[0][2:]
                l_AA = revAA_lookup[left[0:3]]
                l_pos_index = int(left[3:])-1
                right = AA_change_full.split('_')[1].split('ins')[0]
                r_AA = revAA_lookup[right[0:3]]
                r_pos_index = int(right[3:])-1
            
                insertion3 = AA_change_full.split('ins')[1]
                insertion=''
                for i in range(0,len(insertion3),3):
                    insertion=insertion+revAA_lookup[insertion3[i:i+3]]

            elif 'dup' in AA_change_full:
                if '_' not in AA_change_full:
                    left = AA_change_full[2:]
                    l_AA=revAA_lookup[left[0:3]]
                    l_pos_index = int(left[3:-3])
                    r_pos_index = l_pos_index + 1
                    insertion=l_AA

                else:
                    left = AA_change_full.split('_')[0][2:]
                    l_AA = revAA_lookup[left[0:3]]
                    l_pos_index = int(left[3:])-1
                    right = AA_change_full.split('_')[1].split('dup')[0]
                    r_AA = revAA_lookup[right[0:3]]
                    r_pos_index = int(right[3:])-1

                    insertion= ref_dict[t][l_pos_index:r_pos_index+1]
                    #for i in range((r_pos_index-l_pos_index+1)/2):
                    #    insertion=insertion + l_AA + r_AA


            """
            position = re.search(r'[0-9]+',AA_change_full)	
            if position is None:
                print("couldn't find position: {}".format(t))
                print(line)
                bailedOut = True
                continue

            variant_pos = int(position.group())
            """
            if t in blast_dict:

                repeated_substring=False
                one_=False
                two_=False
                variant_has_matched_transcript = True
                mstrgs = blast_dict[t]
                false_blastHits=0
                t_seq = ref_dict[t]
                ref_substr = t_seq[l_pos_index-min(80,l_pos_index):r_pos_index+min(40,len(t_seq)-(r_pos_index))]
                if sum([1 for x in re.finditer(ref_substr, t_seq)]) >1: 
                    unsearchable_transcripts.append((t,AA_change_full,position.group(),ref_substr))
                    
                    continue
                at_least_1_transcript_without_repeats = True
                for mstrg in set(mstrgs):
                    if '_{}:'.format(chrom) not in mstrg: continue                
                    if '1_' in mstrg: one_=True
                    elif '2_' in mstrg: two_=True
                    #t_seq = ref_dict[t]
                    #ref_substr = t_seq[variant_pos-min(7,variant_pos):variant_pos+min(6,len(t_seq)-(variant_pos))]
                    mstrg_MQnovelPeps = []
                    mstrg_seq = proteome_dict[mstrg]
                    alignments=pairwise2.align.localms(ref_substr,mstrg_seq,4,0,-3,-2.5, penalize_end_gaps=(False,True))

                    for a in alignments:
                        align_formatted = pairwise2.format_alignment(*a)
                        print(align_formatted)
                        align_formatted_lst = [x for x in align_formatted.split('\n')]
                        comps = align_formatted_lst[1]
                        ref_align = align_formatted_lst[0]
                        alt_align = align_formatted_lst[2]
                        if sum([1 for x in comps if x=='|']) < 5*(sum([1 for x in comps if x=='.'])+sum([1 for x in comps if x==' '])): 
                            false_blastHits = false_blastHits+1
                            print('false blasthit')
                        
                        gaps = '-'*len(insertion)
                        if gaps not in ref_align: continue
                        ins_index = ref_align.index(gaps)
                        if alt_align[ins_index:ins_index+len(insertion)] == insertion:
                            successfully_recovered_variant=True


                            if mstrg in mstrg_MQnovelPeps_dict:

                                mstrg_pos_start_index = int(alt_align.split(' ')[-2]) - 1
                                mstrg_pos_end_index = mstrg_pos_start_index+len(ref_substr)
                                mstrg_substr = mstrg_seq[mstrg_pos_start_index-min(30,mstrg_pos_start_index):mstrg_pos_end_index+min(30,len(mstrg_seq)-(mstrg_pos_end_index+1))]
                                for pep in mstrg_MQnovelPeps_dict[mstrg]:
                                    pep_align=pairwise2.align.localms(pep,mstrg_seq,4,0,-3,-2.5, penalize_end_gaps=(False,True))[0]
                                    insertion_pos = a[3]+ins_index
                                    pep_start_pos = pep_align[3]
                                    pep_end_pos = pep_align[4]
                                    if pep in mstrg_substr:
                                        mstrg_MQnovelPeps.append(pep)
                                        open(MQnovelPep_mutation_map_outfile,'a').write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pep, AA_change_full, t, enst_uniprot_dict[t], mstrg, mstrg_substr))

                            open(mutation_mstrg_map_outfile,'a').write('{}\t{}\t{}\t{}\t{}\n'.format(AA_change_full, t, enst_uniprot_dict[t], mstrg,alt_align)) 
                            break
                    if len(mstrg_MQnovelPeps)>0: open(mutation_MQevidence_map_outfile,'a').write("{}\t{}\t{}\t{}\t{}\t{}\n".format(AA_change_full, t, enst_uniprot_dict[t], mstrg, alt_align, mstrg_MQnovelPeps))
                    #if successfully_recovered_variant: break


                if successfully_recovered_variant: break
                else: 
                    print(left)
                    print(right) 
                    print(insertion) 
                    if one_ and two_ and len(set(mstrgs)) - false_blastHits>= 2:
                        print('hits: {}; false hits: {}'.format(len(set(mstrgs)), false_blastHits))
                        at_least_1_tr_w_2_isoforms=True
                        true_fail_transcripts.append((t,AA_change_full))
                        true_fail_set.add(t)
                    else: false_fail_transcripts.append((t,AA_change_full))
        if variant_has_matched_transcript: variants_with_expressed_transcripts=variants_with_expressed_transcripts+1
        if successfully_recovered_variant: successes=successes+1
        elif variant_has_matched_transcript and not at_least_1_transcript_without_repeats: unsearchables +=1
        else:
            if at_least_1_tr_w_2_isoforms: true_fails=true_fails+1
            elif variant_has_matched_transcript: fails_1_isoform = fails_1_isoform+1

print('total insertions: {}'.format(total_insertion))
print('insertions with expressed transcripts: {}'.format(variants_with_expressed_transcripts))
print('insertions successfully recovered: {}'.format(successes))
print('fail, with at least 1 transcript with 2 isoforms: {}'.format(true_fails))
print('fail, but without a transcript w 2 isoforms: {}'.format(fails_1_isoform))
print('indeterminable due to repeated sequences in matched transcripts: {}'.format(unsearchables))
print(true_fail_set)
print(true_fail_transcripts)
print(false_fail_transcripts)
print(unsearchable_transcripts)
#print(transcripts_missing_AA)

nonSnpEff_successes=0
for t in transcripts_missing_AA:
	if t in blast_dict:
		successfully_recovered_variant=False
		mstrgs = blast_dict[t]
		for mstrg in mstrgs:
			t_seq = ref_dict[t]
			mstrg_seq = proteome_dict[mstrg]
			alignment=pairwise2.align.localms(t_seq,mstrg_seq,1,-.5,-1,-0.5)[0]
			align_formatted = pairwise2.format_alignment(*alignment)
			print(align_formatted)
			align_formatted_lst = [x[2:] for x in align_formatted.split('\n')]
			comps = align_formatted_lst[1]
			mismatch_indices = [i for i,val in enumerate(comps) if val == '.']
			ref_align = align_formatted_lst[0]
			alt_align = align_formatted_lst[2]
			if len(comps)>40 and sum([1 for x in comps if x=='.'])>=1 and sum([1 for x in comps if x=='|'])/sum([1 for x in comps if x=='.'])>20:
				for x in comps:
					if x=='.': nonSnpEff_successes=nonSnpEff_successes+1
				successfully_recovered_variant=True
				print('success!')
			"""
			for i in mismatch_indices:
				if ref_align[i] == aa_orig and alt_align[i] == aa_new:
					successfully_recovered_variant=True
					print('success!')
					break
			"""
			if successfully_recovered_variant: break
		#if successfully_recovered_variant: nonSnpEff_successes=nonSnpEff_successes+1

print('non-snpEff successes: {}'.format(nonSnpEff_successes))

