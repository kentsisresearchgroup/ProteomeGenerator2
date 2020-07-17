import sys, re, difflib
from Bio import pairwise2

proteome_infile = sys.argv[1]
vcf_infile = sys.argv[2]
blast_outfmt6_infile = sys.argv[3]
ref_fasta_file = sys.argv[4]

novelPep_mstrgs_infile = sys.argv[5]

mutation_mstrg_map_outfile = sys.argv[6]
mutation_MQevidence_map_outfile = sys.argv[7]
MQnovelPep_mutation_map_outfile = sys.argv[8]

blast_dict=dict()
enst_uniprot_dict=dict()
ref_blastHits=set()
blast_mismatch_rate_dict=dict()
with open(blast_outfmt6_infile) as blast:
        for hit in blast:
                entry = hit.split('\t')
                mstrg = entry[0][3:-1]
                ref_header = entry[1]
                enst = re.search(r'ENST[0-9]+',ref_header).group()
                enst_uniprot_dict[enst] = ref_header.split('|')[5]
                mismatch_rate = float(int(entry[4])/int(entry[3]))
                if mismatch_rate > 0.05: continue
                if enst in blast_dict:
                        lst = blast_dict[enst]
                        lst.append(mstrg)
                        blast_dict[enst] = lst
                else:
                        blast_dict[enst] = [mstrg]
                blast_mismatch_rate_dict[(enst,mstrg)] = mismatch_rate


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

frameshift_dict=dict()
with open(vcf_infile,'r') as f:
    for line in f:
        match = re.search(r'frameshift',line)
        if match is None: continue
        transcripts = re.findall(r'ENST[0-9]+',line)
        for t in transcripts: 
            if t not in frameshift_dict: frameshift_dict[t]=1
            else: frameshift_dict[t] += 1
        

total_frameshift = 0
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

with open(vcf_infile,'r') as f:
    for line in f:
        chrom = line.split('\t')[0]
        successfully_recovered_variant = False
        at_least_1_tr_w_2_isoforms=False
        at_least_1_transcript_without_repeats=False
        variant_has_matched_transcript = False

        match = re.search(r'frameshift',line)
        if match is None: continue
        total_frameshift = total_frameshift+1

        # iterate over all predicted coding transcripts in the SnpEff variant entry
        transcripts = re.findall(r'ENST[0-9]+',line)
        for t in transcripts:
            
            # proceed only if coding transcript
            in_cds = re.search(r'frameshift',line.split(t)[0].split(',')[-1])
            if not in_cds: continue
            if t not in ref_dict:
                print('{}: transcript not in ref protein fasta'.format(t))
                continue
            print(t)

            # extract amino acid change prediction
            AA_change_full = re.search(r'p\.[A-Z][a-z][a-z][0-9]+fs',line.split(t)[1].split(',')[0])
            if AA_change_full is None:
                print("couldn't find AA change: {}".format(t))
                transcripts_missing_AA.add(t)
                bailedOut = True
                continue
            AA_change_full = AA_change_full.group()
            print(AA_change_full)

            aa_new = revAA_lookup[re.findall(r'[A-Z][a-z][a-z]',AA_change_full)[0]]

            # extract positional change prediction
            position = int(re.findall(r'[0-9]+',AA_change_full)[0])
            if position is None:
                print("couldn't find position: {}".format(t))
                print(line)
                bailedOut = True
                continue
            variant_pos_index = position-1

            if t in blast_dict: # proceed if the predicted ref transcript is present in the proteome
                variant_has_matched_transcript = True
                mstrgs = blast_dict[t] # find proteome isoforms matching the ref transcript per blast

                t_seq = ref_dict[t]
                ref_substr = t_seq[variant_pos_index-min(25,variant_pos_index):variant_pos_index+min(15,len(t_seq)-(variant_pos_index+1))]

                one_=False # "failed" variant recoveries, in transcripts for which there are not copies in both haplotype 1 and haplotype 2, will not be counted due to the possibility of the mutant copy being silenced, mutant exon being skipped, etc
                two_=False
                false_blastHits=0 # not all "hits" from blast are valid
                
                # "failed" variant recoveries, in transcripts  where the reference substring of interest appears more than once, will not be counted due to the possibility of aligning to the wrong instance of the repeated substirng
                if sum([1 for x in re.finditer(ref_substr, t_seq)]) >1: 
                    unsearchable_transcripts.append((t,AA_change_full,position,ref_substr))
                    continue
                at_least_1_transcript_without_repeats = True

                for mstrg in set(mstrgs):
                    if '_{}:'.format(chrom) not in mstrg: continue

                    mstrg_seq = proteome_dict[mstrg]
                    mstrg_MQnovelPeps = []    

                    alignments=pairwise2.align.localms(ref_substr,mstrg_seq,4,0,-3,-2.5, penalize_end_gaps=(False,True))

                    for a in alignments:
                        align_formatted = pairwise2.format_alignment(*a)
                        print(align_formatted)
                        align_formatted_lst = [x for x in align_formatted.split('\n')]
                        ref_align = align_formatted_lst[0].split(' ')[-1]
                        alt_align = align_formatted_lst[2].split(' ')[-1]
                        alt_align_start_pos = int(re.search(r'[0-9]+',align_formatted_lst[2]).group())
                        comps = align_formatted_lst[1][-1*len(ref_align):]
                        mismatch_indices = [i for i,val in enumerate(comps) if val == '.']
                        """if sum([1 for x in comps if x=='|']) < 3*(sum([1 for x in comps if x=='.'])+sum([1 for x in comps if x==' '])): 
                            false_blastHits = false_blastHits+1
                            print('false blasthit') """
                        if '1_' in mstrg: one_=True
                        elif '2_' in mstrg: two_=True
                        for i in mismatch_indices:
                            print('{}\t{}\t{}'.format(i,ref_align[i],alt_align[i]))
                            pre_i_matches = sum([1 for x in comps[0:i] if x == '|'])
                            pre_i_mismatches = sum([1 for x in comps[0:i] if x == '.' or x == ' '])
                            post_i_matches = sum([1 for x in comps[i:] if x == '|'])
                            post_i_mismatches = sum([1 for x in comps[i:] if x == '.' or x == ' '])
                            #if alt_align[i] == aa_new: print('{}\t{}\t{}\t{}'.format(i,new_stop-1,alt_align_start_pos+i+new_stop-1,len(mstrg_seq)))
                            if (alt_align[i] == aa_new and ((pre_i_matches>5*pre_i_mismatches and post_i_mismatches>post_i_matches) or (pre_i_mismatches>pre_i_matches and post_i_matches > 3*post_i_mismatches))) or (len(ref_substr) and pre_i_matches>5*pre_i_mismatches and alt_align_start_pos+(len(ref_substr)-1)<len(mstrg_seq)):
                            #if ref_align[i] == aa_orig and alt_align[i] == aa_new and (abs(len(mstrg_seq) - (alt_align_start_pos+i+(new_stop-1)))<=1 or frameshift_dict[t]>1):
                                successfully_recovered_variant=True
                                
                                if mstrg in mstrg_MQnovelPeps_dict:
                                
                                    #mstrg_pos_start_index = int(alt_align.split(' ')[0]) - 1
                                    mstrg_pos_start_index = a[3]
                                    print(mstrg_pos_start_index)
                                    mstrg_pos_end_index = mstrg_pos_start_index+len(ref_substr) 
                                    mstrg_substr = mstrg_seq[mstrg_pos_start_index-min(20,mstrg_pos_start_index):mstrg_pos_end_index+min(20,len(mstrg_seq)-(mstrg_pos_end_index+1))]
                                    for pep in mstrg_MQnovelPeps_dict[mstrg]:
                                        if pep in mstrg_substr:
                                            mstrg_MQnovelPeps.append(pep)
                                            open(MQnovelPep_mutation_map_outfile,'a').write("{}\t{}\t{}\t{}\t{}\t{}\n".format(pep, AA_change_full, t, enst_uniprot_dict[t], mstrg, mstrg_substr))
                                        
                                open(mutation_mstrg_map_outfile,'a').write('{}\t{}\t{}\t{}\t{}\n'.format(AA_change_full, t, enst_uniprot_dict[t], mstrg,alt_align))
                                #open(mutation_mstrg_map_outfile,'a').write('{}\t{}\t{}\n'.format(AA_change_full,mstrg.split('|')[1],alt_align))
                                break
                        if successfully_recovered_variant: break
                    if len(mstrg_MQnovelPeps)>0: open(mutation_MQevidence_map_outfile,'a').write("{}\t{}\t{}\t{}\t{}\t{}\n".format(AA_change_full, t, enst_uniprot_dict[t], mstrg, alt_align, mstrg_MQnovelPeps))
                    #if successfully_recovered_variant: break


                if successfully_recovered_variant: break
                else: 
                    #print(aa_orig)
                    print(aa_new)
                    print(position) 
                    if one_ and two_ and len(set(mstrgs)) - false_blastHits>= 2:
                        print('hits: {}; false hits: {}'.format(len(mstrgs), false_blastHits))
                        at_least_1_tr_w_2_isoforms=True
                        true_fail_transcripts.append((t,AA_change_full,position))
                        true_fail_set.add(t)
                    else: false_fail_transcripts.append((t,AA_change_full,position))
        if variant_has_matched_transcript: variants_with_expressed_transcripts=variants_with_expressed_transcripts+1
        if successfully_recovered_variant: successes=successes+1
        elif variant_has_matched_transcript and not at_least_1_transcript_without_repeats: unsearchables +=1
        else:
            if at_least_1_tr_w_2_isoforms: true_fails=true_fails+1
            elif variant_has_matched_transcript: fails_1_isoform = fails_1_isoform+1

print('total frameshifts: {}'.format(total_frameshift))
print('frameshifts with expressed transcripts: {}'.format(variants_with_expressed_transcripts))
print('frameshifts successfully recovered: {}'.format(successes))
print('fail, with at least 1 transcript with 2 isoforms: {}'.format(true_fails))
print('fail, but without a transcript w 2 isoforms: {}'.format(fails_1_isoform))
print('indeterminable due to repeated sequences in matched transcripts: {}'.format(unsearchables))
print(true_fail_set)
print(true_fail_transcripts)
print(false_fail_transcripts)
print(unsearchable_transcripts)
#print(transcripts_missing_AA)

