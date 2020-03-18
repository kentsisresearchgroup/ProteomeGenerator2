import sys
#fasta_path = sys.argv[1]
out_file = sys.argv[1]
fasta = {}
for filename in sys.argv[2:]:
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                #if active_sequence_name not in fasta:
                #     fasta[active_sequence_name] = ''
                if active_sequence_name in fasta:
                    for i in range(2,500):
                        if '{}_{}'.format(active_sequence_name,i) not in fasta:
                            active_sequence_name = '{}_{}'.format(active_sequence_name,i)
                            break
                fasta[active_sequence_name] = ''
                continue
            sequence = line
            fasta[active_sequence_name] = fasta[active_sequence_name] + sequence
for k,v in fasta.items():
	open('check_fasta.txt','a').write('{}\n{}\n'.format(k,v))
reverse_dict = {}
for header, seq in fasta.items():
    if seq not in reverse_dict:
        reverse_dict[seq] = header    
    else:
        if 'ENST' not in header and 'MSTRG' not in header:
            reverse_dict[seq] = header
        elif 'ENST' in header:
            reverse_dict[seq] = header

import textwrap
with open(out_file,'w') as out:
    for seq,header in reverse_dict.items():
        out.write('>pg|{}|\n'.format(header))
        #out.write('{}\n'.format(seq))
        #[out.write('{}\n'.format(x)) for x in textwrap.wrap(seq, 80)]
        out.write('{}\n'.format('\n'.join(textwrap.wrap(seq, 80))))
#fasta
