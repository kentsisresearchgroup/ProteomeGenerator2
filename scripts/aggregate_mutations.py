import sys,re,os

pathname=os.path.dirname(sys.argv[1])
mutation_type = re.search(r'(chr[0-9X]+).([a-z]+).([a-z]+)',os.path.basename(sys.argv[1])).group(2)
mutation_maps = [x for x in sys.argv if '.{}.map'.format(mutation_type) in x]
novelpep_maps = [x for x in sys.argv if '.novelPep_{}.map'.format(mutation_type) in x]
mqevidence_maps = [x for x in sys.argv if 'MQevidence' in x]

mutation_mstrg_dict=dict()
novelpep_mutation_dict=dict()
mutation_novelpep_dict=dict()

for _map in mutation_maps:
	with open(_map) as f:
		for line in f:
			entry = line.strip().split('\t')
			key_tuple = ((entry[1],entry[0]),entry[2])
			mstrg_fmt=re.sub(r'([0-9]_chr)',r' \1',entry[3])
			if key_tuple not in mutation_mstrg_dict:
				mutation_mstrg_dict[key_tuple] = [mstrg_fmt]
			else: 
				mstrgs=mutation_mstrg_dict[key_tuple]
				mstrgs.append(mstrg_fmt)
				mutation_mstrg_dict[key_tuple] = mstrgs

for _map in novelpep_maps:
	with open(_map) as f:
		for line in f:
			entry = line.strip().split('\t')			 
			key_tuple = (entry[0],(entry[1],entry[2]),entry[3])
			mstrg_fmt=re.sub(r'([0-9]_chr)',r' \1',entry[4])
			if key_tuple not in novelpep_mutation_dict:
                                novelpep_mutation_dict[key_tuple] = [mstrg_fmt]
			else:
                                mstrgs=novelpep_mutation_dict[key_tuple]
                                mstrgs.append(mstrg_fmt)
                                novelpep_mutation_dict[key_tuple] = mstrgs

for _map in mqevidence_maps:
	with open(_map) as f:
		for line in f:
			entry = line.strip().split('\t')
			key_tuple = ((entry[1],entry[0]),entry[2])
			mstrg_fmt=re.sub(r'([0-9]_chr)',r' \1',entry[3])
			if key_tuple not in mutation_novelpep_dict:
				mutation_novelpep_dict[key_tuple] = [mstrg_fmt,entry[5]]

with open(os.path.join(pathname,'combined.{}.map'.format(mutation_type)),'w') as f:
	for k in mutation_mstrg_dict:
		(transcript_tup,prot) = k
		mstrgs = mutation_mstrg_dict[k]
		f.write('{}\t{}\t{}\n'.format(transcript_tup,prot,mstrgs))

with open(os.path.join(pathname,'combined.novelPep_{}.map'.format(mutation_type)),'w') as f:
	for k in novelpep_mutation_dict:
		(peptide,transcript_tup,prot) = k
		mstrgs = novelpep_mutation_dict[k]
		f.write('{}\t{}\t{}\t{}\n'.format(peptide,transcript_tup,prot,mstrgs))

with open(os.path.join(pathname,'combined.{}_MQevidence.map'.format(mutation_type)),'w') as f:
	for k in mutation_novelpep_dict:
		(transcript_tup,prot) = k
		v = mutation_novelpep_dict[k]
		f.write('{}\t{}\t{}\t{}\n'.format(transcript_tup,prot,v[1],v[0]))
