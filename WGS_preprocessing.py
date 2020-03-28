import os
PG2_HOME=config['directories']['PG2_installation_dir']
WD=config['directories']['working_and_output_dir']
workdir: WD
TMP=config['directories']['optional']['temp_dir']

STOCK_GENOME_GTF=config['stock_references']['genome']['gtf']
STOCK_GENOME_FASTA=config['stock_references']['genome']['fasta']
prebuilt_bwa = (config['stock_references']['genome']['optional_aligner_indices']['BWA']['BWA_index_dir'] is not None)
if prebuilt_bwa:
    BWA_INDEX_DIR=config['stock_references']['genome']['optional_aligner_indices']['BWA']['BWA_index_dir']
    BWA_PREFIX=config['stock_references']['genome']['optional_aligner_indices']['BWA']['BWA_prefix']
    BWA_INDEX=os.path.join(BWA_INDEX_DIR,BWA_PREFIX+'.sa')
    STOCK_GENOME_FASTA = os.path.join(BWA_INDEX_DIR,BWA_PREFIX.strip('.fa')+'.fa')
else:
    BWA_INDEX_DIR=os.path.dirname(STOCK_GENOME_FASTA)
    BWA_INDEX=STOCK_GENOME_FASTA+'.sa'


snakemake.utils.makedirs('out/benchmarks')
snakemake.utils.makedirs('out/logs')
snakemake.utils.makedirs('out/readgroups')

input_file_format = config['input_files']['genome_personalization_module']['input_file_format']


if input_file_format == 'bam':
    SAMPLES=config['input_files']['genome_personalization_module']['bam_inputs'].keys()
    """
    rule wgs_00bam_PrepareRGwiseOutputMapForBamRevert:
        output: tsv="out/WGS/{tumor_or_normal}/{sample}.reverted_ubam_RGmap.tsv"
        params: n="1", R="'span[hosts=1] rusage[mem=4]'", o="out/logs/rgmap.out", eo="out/logs/rgmap.err", J="create_rgmap", RG_header_list= lambda wildcards: config['input_files']['genome_personalization_module']['bam_inputs'][wildcards.sample]['read_groups']
        run:
            with open(output.tsv,'w') as f:
                f.write('READ_GROUP_ID\tOUTPUT\n')
                for header in params.RG_header_list:
                    rg_id = [x for x in header.split(' ') if 'ID:' in x][0].replace('ID:','')
                    f.write('{}\t{}/out/WGS/{}/{}.{}.unmapped.bam\n'.format(rg_id,WD,wildcards.tumor_or_normal,wildcards.sample,rg_id))
    """

    rule wgs_01bam_RevertToUnmappedBAM:
        input: bam=lambda wildcards: config['input_files']['genome_personalization_module']['bam_inputs'][wildcards.sample]['bam_file'], tsv="out/WGS/{tumor_or_normal}/{sample}.reverted_ubam_RGmap.tsv",ref=STOCK_GENOME_FASTA
        output: temp("out/WGS/{tumor_or_normal}/{sample}.RG.unmapped.bam")
        params: n="16", R="'span[hosts=1] rusage[mem=10]'", o="out/logs/revert_bam.out", eo="out/logs/revert_bam.err", J="revert_bam"
        conda: "envs/bwa_picard_samtools.yaml"
        shell: "picard -Xmx160g RevertSam INPUT={input.bam} TMP_DIR={TMP} R={input.ref} \
                  MAX_RECORDS_IN_RAM=25000000 \
                  O={output} \
                  MAX_DISCARD_FRACTION=0.03 VALIDATION_STRINGENCY=SILENT"

elif input_file_format == 'fastq':

    """   
    rule wgs_00_CleanRawFastqsWithFastp:
        input: read_one=lambda wildcards: config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R1_fq.gz'],read_two=lambda wildcards: config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R2_fq.gz']
        output: read_one=temp("out/WGS/{tumor_or_normal}/{sample}.{readgroup}.fastp.1.fq"), read_two=temp("out/WGS/{tumor_or_normal}/{sample}.{readgroup}.fastp.2.fq")
        params: n="8", R="'span[hosts=1] rusage[mem=4]'", o="out/logs/fastp.out", eo="out/logs/fastp.err", J="fastp_qc"
        conda: "envs/bwa_picard_samtools.yaml"
        shell: "fastp -i {input.read_one} -o {output.read_one} -I {input.read_two} -O {output.read_two}"
    """   
    rule wgs_01fq_GenerateRGWiseUbamsFromFq:
        #input: read_one="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.fastp.1.fq", read_two="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.fastp.2.fq"
        input: read_one=lambda wildcards: config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R1_fq.gz'],read_two=lambda wildcards: config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R2_fq.gz']
        output: temp("out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.bam")
        params: n="8", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/fastq2ubam.out", eo="out/logs/fastq2ubam.err", J="fq2ubams", \
                tumor_or_normal=lambda wildcards: config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['matched_sample_params']['tumor_or_normal'] if config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['matched_sample_params']['is_matched_sample'] else 'tumor', \
                output_bamfile=lambda wildcards: WD + "/out/WGS/{}/{}.{}.unmapped.bam".format(wildcards.tumor_or_normal,wildcards.sample,wildcards.readgroup), \
                library_name = lambda wildcards: config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['optional']['library_id'], sequencing_platform= lambda wildcards:config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['optional']['sequencing_platform'], \
                max_records_in_ram = 250000*64
        conda: "envs/bwa_picard_samtools.yaml"
        shell:  "picard FastqToSam -Xmx96g TMP_DIR={TMP} \
                   FASTQ={input.read_one} FASTQ2={input.read_two} OUTPUT={params.output_bamfile} \
                   READ_GROUP_NAME={wildcards.readgroup} SAMPLE_NAME={wildcards.sample} \
                   LIBRARY_NAME={params.library_name} PLATFORM={params.sequencing_platform} \
                   MAX_RECORDS_IN_RAM={params.max_records_in_ram} SORT_ORDER=queryname"
    """
    rule wgs_02_NameSortedUbam2Fastq4Alignment:
        input: ubam_file="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.bam"
        output: fq1=temp("out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.1.fq"), fq2=temp("out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.2.fq")
        params: n="4",  R="'span[hosts=1]'", o="out/logs/bam2fq.out", eo="out/logs/bam2fq.err", J="bam2fq"
        conda: "envs/bwa_picard_samtools.yaml"
        shell: "picard SamToFastq I={input.ubam_file} FASTQ={output.fq1} SECOND_END_FASTQ={output.fq2}"

    rule wgs_03_BwaAndSortAndMergeBamAlignment_FQ:
        input:  ubam="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.bam",read_one="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.1.fq", read_two="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.2.fq", ref_idx=STOCK_GENOME_FASTA+".fai", ref_dict=STOCK_GENOME_FASTA.strip('fa')+'dict',bwa_idx=BWA_INDEX
        output: "out/WGS/{tumor_or_normal}/{sample}.{readgroup}.aligned_sorted_ubam-merged.bam"
        params: n="32", R="'span[hosts=1] rusage[mem=2]'", o="out/logs/bwa_n_mergebams.out", eo="out/logs/bwa_n_mergebams.err", J="bwa_mergebams", tumor_or_normal=lambda wildcards: config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['matched_sample_params']['tumor_or_normal'] if config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['matched_sample_params']['is_matched_sample'] else 'tumor', output_bamfile=lambda wildcards:WD + "/out/WGS/{}/{}.{}.aligned_sorted.bam".format(wildcards.tumor_or_normal,wildcards.sample,wildcards.readgroup)
        conda: "envs/bwa_picard_samtools.yaml"
        shell: "bwa mem -M -t {params.n} {STOCK_GENOME_FASTA} {input.read_one} {input.read_two} | samtools view -Shu -@ {params.n} - | \
                picard MergeBamAlignment -Xmx64g ALIGNED_BAM=/dev/stdin UNMAPPED_BAM={input.ubam} \
                  REFERENCE_SEQUENCE={STOCK_GENOME_FASTA} OUTPUT={output} CREATE_INDEX=true \
                  CLIP_OVERLAPPING_READS=true \
                  INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
                  ATTRIBUTES_TO_RETAIN=XS \
                  MAX_RECORDS_IN_RAM=30000000 TMP_DIR={TMP} \
                  SORT_ORDER=queryname"
    """

rule wgs_02_NameSortedUbam2Fastq4Alignment:
    input: ubam_file="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.bam"
    output: fq1=temp("out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.1.fq.gz"), fq2=temp("out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.2.fq.gz")
    params: n="4",  R="'span[hosts=1]'", o="out/logs/bam2fq.out", eo="out/logs/bam2fq.err", J="bam2fq", \
            tmp_fq1=os.path.join(TMP,"{sample}.{readgroup}.unmapped.ubam2fq.1.fq"),tmp_fq2=os.path.join(TMP,"{sample}.{readgroup}.unmapped.ubam2bq.2.fq")
    conda: "envs/bwa_picard_samtools.yaml"
    shell: "picard SamToFastq I={input.ubam_file} FASTQ={params.tmp_fq1} SECOND_END_FASTQ={params.tmp_fq2}; gzip -c {params.tmp_fq1} > {output.fq1}; gzip -c {params.tmp_fq2} > {output.fq2}"

rule wgs_03_BwaAndSortAndMergeBamAlignment_FQ:
    input:  ubam="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.bam",read_one="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.1.fq.gz", read_two="out/WGS/{tumor_or_normal}/{sample}.{readgroup}.unmapped.ubam2fq.2.fq.gz", ref_idx=STOCK_GENOME_FASTA+".fai", ref_dict=os.path.splitext(STOCK_GENOME_FASTA)[0]+'.dict',bwa_idx=BWA_INDEX
    output: "out/WGS/{tumor_or_normal}/{sample}.{readgroup}.aligned_sorted_ubam-merged.bam"
    params: n="32", R="'span[hosts=1] rusage[mem=2]'", o="out/logs/bwa_n_mergebams_{sample}_{readgroup}.out", eo="out/logs/bwa_n_mergebams_{sample}_{readgroup}.err", J="bwa_mergebams", tumor_or_normal=lambda wildcards: config['input_files']['genome_personalization_module'][input_file_format+'_inputs'][wildcards.sample]['matched_sample_params']['tumor_or_normal'] if config['input_files']['genome_personalization_module'][input_file_format+'_inputs'][wildcards.sample]['matched_sample_params']['is_matched_sample'] else 'tumor', output_bamfile=lambda wildcards:WD + "/out/WGS/{}/{}.{}.aligned_sorted.bam".format(wildcards.tumor_or_normal,wildcards.sample,wildcards.readgroup)
    conda: "envs/bwa_picard_samtools.yaml"
    shell: "bwa mem -M -t {params.n} {STOCK_GENOME_FASTA} {input.read_one} {input.read_two} | samtools view -Shu -@ {params.n} - | \
            picard MergeBamAlignment -Xmx64g ALIGNED_BAM=/dev/stdin UNMAPPED_BAM={input.ubam} \
              REFERENCE_SEQUENCE={STOCK_GENOME_FASTA} OUTPUT={output} CREATE_INDEX=true \
              CLIP_OVERLAPPING_READS=true \
              INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
              ATTRIBUTES_TO_RETAIN=XS \
              MAX_RECORDS_IN_RAM=30000000 TMP_DIR={TMP} \
              SORT_ORDER=queryname"

if input_file_format == 'fastq':
    rule wgs_04_MergeAllRGs:
        input: lambda wildcards: expand("out/WGS/{{tumor_or_normal}}/{{sample}}.{readgroup}.aligned_sorted_ubam-merged.bam", readgroup=config['input_files']['genome_personalization_module']['fastq_inputs'][wildcards.sample]['read_groups'].keys())
        output: temp("out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged.bam")
        params: n="16", R="'span[hosts=1] rusage[mem=4]'", o="out/logs/merge_sample_RGs.out", eo="out/logs/merge_sample_RGs.err", J="merge_sample_RGs"
        conda: "envs/bwa_picard_samtools.yaml"
        shell: "picard MergeSamFiles -Xmx64g \
                  $(echo '{input}' | sed -r \'s/[^ ]+/INPUT=&/g') \
                  OUTPUT={output} SORT_ORDER=queryname TMP_DIR={TMP}"

elif input_file_format == 'bam':
    rule wgs_04_MergeAllRGs:
        #input: lambda wildcards: expand("out/WGS/{{tumor_or_normal}}/{{sample}}.{readgroup}.aligned_sorted_ubam-merged.bam", readgroup=[[q for q in x.split(' ') if 'ID:' in q][0].replace('ID:','') for x in config['input_files']['genome_personalization_module']['bam_inputs'][wildcards.sample]['read_groups']])
        input: "out/WGS/{tumor_or_normal}/{sample}.RG.aligned_sorted_ubam-merged.bam"
        output: temp("out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged.bam")
        params: n="16", R="'span[hosts=1] rusage[mem=4]'", o="out/logs/merge_sample_RGs.out", eo="out/logs/merge_sample_RGs.err", J="merge_sample_RGs"
        conda: "envs/bwa_picard_samtools.yaml"
        shell: "picard MergeSamFiles -Xmx64g \
                  $(echo '{input}' | sed -r \'s/[^ ]+/INPUT=&/g') \
                  OUTPUT={output} SORT_ORDER=queryname TMP_DIR={TMP}"

rule wgs_05_MarkDuplicatesAndCoordSort:
    input: "out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged.bam"
    output: bam=temp("out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup.bam"), metrics="out/WGS/{tumor_or_normal}/{sample}.dup_metrics.txt" # TODO: remove 'mates-fixed' from this and subsequent input/outputs
    params: n="24", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/markdups.out", eo="out/logs/markdups.err", J="mark_dups"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx32g' MarkDuplicates \
              -I {input} -ASO queryname -O /dev/stdout -M {output.metrics} \
              --TMP_DIR {TMP} --MAX_RECORDS_IN_RAM 8000000 | \
            gatk --java-options '-Xmx155g' SortSam \
              -I /dev/stdin -O {output.bam} --SORT_ORDER coordinate \
              --TMP_DIR {TMP} --MAX_RECORDS_IN_RAM 30000000"

rule wgs_06_IndexBam:
    input: "out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup.bam"
    output: "out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup.bai"
    params: n="4", R="'span[hosts=1] rusage[mem=4]'", o="out/logs/index_bam.out", eo="out/logs/index_bam.err", J="index_bam"
    conda: "envs/bwa_picard_samtools.yaml"
    shell: "samtools index {input} {output}"

rule wgs_00a_CreateIntervalsWithoutNs:
    input: STOCK_GENOME_FASTA
    output: "out/WGS/no_N.interval_list"
    params: n="1", R="'span[hosts=1] rusage[mem=6]'", \
            o="out/logs/create_intervals.out", eo="out/logs/create_intervals.err", \
            J="create_intervals"
    conda: "envs/gatk4.yaml"
    shell: "picard ScatterIntervalsByNs R={input} O={output}"

NUM_BQSR_INTERVALS=config['parameters']['genome_personalization_module']['pre-processing']['advanced']['BQSR_intervals_scatter']
rule wgs_00b_ScatterBaseRecalibratorIntervals:
    input: "out/WGS/no_N.interval_list"
    output: expand("out/WGS/intervals/BQSR/{interval}-scattered.interval_list",interval=[str(x).zfill(4) for x in range(NUM_BQSR_INTERVALS)])
    params: n="1", R="'span[hosts=1] rusage[mem=4]'", \
            o="out/logs/intervals/BQSR_scatter.out", eo="out/logs/intervals/BQSR_scatter.err", \
            J="scatter_intervals_BQSR"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx4g' SplitIntervals \
              -R {STOCK_GENOME_FASTA} -L {input} -scatter {NUM_BQSR_INTERVALS} \
              -O out/WGS/intervals/BQSR"


snakemake.utils.makedirs('out/BQSR')
snakemake.utils.makedirs('out/logs/BQSR')

snp_known_sites=config['parameters']['genome_personalization_module']['variant_calling']['resources']['germline']['snps_db']
indel_known_sites=config['parameters']['genome_personalization_module']['variant_calling']['resources']['germline']['indels_db']
rule wgs_07_BaseRecalibrator:
    input: bam="out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup.bam", bai="out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup.bai",interval_list="out/WGS/intervals/BQSR/{interval}-scattered.interval_list"
    output: temp("out/WGS/{tumor_or_normal}/BQSR/scattered-calculate/{sample}.bqsr-{interval}.report")
    params: n="8", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/BQSR/{interval}_BQSR.out", eo="out/logs/BQSR/{interval}_BQSR.err", \
	    J="base_BQSR_{interval}"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx64g' BaseRecalibrator -I {input.bam} -R {STOCK_GENOME_FASTA} --known-sites {snp_known_sites} --known-sites {indel_known_sites} -O {output} -L {input.interval_list}"

rule wgs_08_GatherBqsrReports:
    input: expand("out/WGS/{{tumor_or_normal}}/BQSR/scattered-calculate/{{sample}}.bqsr-{interval}.report", interval=[str(x).zfill(4) for x in range(NUM_BQSR_INTERVALS)])
    output: "out/WGS/{tumor_or_normal}/BQSR/{sample}.bqsr-calculated.report"
    params: n="1", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/gather_reports.out", eo="out/logs/gather_reports.err", \
	    J="gather_reports"
    conda: "envs/gatk4.yaml"
    shell: "gatk GatherBQSRReports $(echo {input} | sed -r 's/[^ ]+/-I &/g') -O {output}"

rule wgs_09a_ApplyBQSRToMappedReads:
    input: bqsr="out/WGS/{tumor_or_normal}/BQSR/{sample}.bqsr-calculated.report",bam="out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup.bam",interval_list="out/WGS/intervals/BQSR/{interval}-scattered.interval_list"
    output: temp("out/WGS/{tumor_or_normal}/BQSR/scattered-apply/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_{interval}_scatteredBQSR.bam")
    wildcard_constraints: interval="\d+"
    params: n="2", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/mapped_bqsr_{interval}.out", eo="out/logs/mapped_bqsr{interval}.err", \
	    J="mapped_bqsr"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx8g' ApplyBQSR \
	      -R {STOCK_GENOME_FASTA} -I {input.bam} -O {output} \
	      --bqsr-recal-file {input.bqsr} \
	      -L {input.interval_list}"

rule wgs_09b_ApplyBQSRToUnmappedReads:
    input: bqsr="out/WGS/{tumor_or_normal}/BQSR/{sample}.bqsr-calculated.report",bam="out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup.bam"
    output: temp("out/WGS/{tumor_or_normal}/BQSR/scattered-apply/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_unmapped_scatteredBQSR.bam")
    params: n="2", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/unmapped_bqsr.out", eo="out/logs/unmapped_bqsr.err", \
	    J="unmapped_bqsr"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx8g' ApplyBQSR \
	      -R {STOCK_GENOME_FASTA} -I {input.bam} -O {output} \
	      --bqsr-recal-file {input.bqsr} \
	      -L unmapped"
	     
rule wgs_10_GatherRecalibratedBAMs:
    input: mapped=expand("out/WGS/{{tumor_or_normal}}/BQSR/scattered-apply/{{sample}}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_{interval}_scatteredBQSR.bam",interval=[str(x).zfill(4) for x in range(NUM_BQSR_INTERVALS)]), \
	   unmapped="out/WGS/{tumor_or_normal}/BQSR/scattered-apply/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_unmapped_scatteredBQSR.bam"
    output: bam="out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bam",bai="out/WGS/{tumor_or_normal}/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bai"
    params: n="4", R="'span[hosts=1] rusage[mem=8]'", \
	    o="out/logs/gather_recalibrated_bams.out", eo="out/logs/gather_recalibrated_bams.err", \
	    J="gather_recalibrated_bams"
    conda: "envs/gatk4.yaml"
    shell: "picard MergeSamFiles -Xmx32g \
              $(echo '{input.mapped}' | sed -r \'s/[^ ]+/INPUT=&/g') \
	      INPUT={input.unmapped} OUTPUT={output.bam} \
	      CREATE_INDEX=true CREATE_MD5_FILE=true"

rule util_CreateGenomeDict:
    input: STOCK_GENOME_FASTA
    output: os.path.splitext(STOCK_GENOME_FASTA)[0]+'.dict'
    params: n="1", R="'span[hosts=1] rusage[mem=18]'", o="out/logs/create_refDict.out", eo="out/logs/create_refDict.err", J="create_refDict"
    conda: "envs/gatk4.yaml"
    shell: "picard -Xmx16g CreateSequenceDictionary R={input} O={output}"

rule util_CreateGenomeFastaIndex:
    input: STOCK_GENOME_FASTA
    output: STOCK_GENOME_FASTA+".fai"
    params: n="1", R="'span[hosts=1] rusage[mem=18]'", o="out/logs/create_refIdx.out", eo="out/logs/create_refDict.err", J="create_refIdx"
    conda: "envs/gatk4.yaml"
    shell: "samtools faidx {input}"

rule util_CreateBWAindex:
    input: STOCK_GENOME_FASTA
    output: BWA_INDEX
    params: n="24", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/bwa_index.out", eo="out/logs/bwa_index.err", J="create_bwa_index"
    conda: "envs/bwa_picard_samtools.yaml"
    shell: "bwa index -a bwtsw {input}"
