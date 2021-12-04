shell.executable("/bin/bash")

import os, multiprocessing
import os.path

### Begin User Variables ###

#Initialize home/working/output directories
PG2_HOME = config['directories']['PG2_installation_dir']
WD = config["directories"]["working_and_output_dir"]
workdir: WD
TMP = config["directories"]['optional']["temp_dir"]

# Initialize stock references
STOCK_GENOME_FASTA = config['stock_references']['genome']['fasta']
STOCK_GENOME_GTF = config['stock_references']['genome']['gtf']
STOCK_PROTEOME_FASTA = config['stock_references']['proteome']['fasta']

# Cohort/Organism info
COHORT = config['input_files']['genome_personalization_module']['cohort_or_organism_name']

# User-defined parameters
ORF = config['parameters']['protein_prediction_main']['min_peptide_ORF_length']
nuc_ORF=str(int(ORF)*3)

### End User Variables ###


### Workflow Control ###

# Boolean workflow switches
creating_custom_genome = config['user_defined_workflow']['genome_personalization_module']['enable_genome_personalization_module'] # boolean to check if custom reference is to be created
RNA_seq_module_enabled = config['user_defined_workflow']['RNA-seq_module']['enable_RNA-seq_module']

continuing_after_genome_personalization = config['user_defined_workflow']['genome_personalization_module']['continuation']['just_created_custom_genome']
assert not (creating_custom_genome and continuing_after_genome_personalization), "The genome personalization module and 'just created custom genome' flag (i.e. already ran genome personalization) are both true. Please disable at least one in the configuration file"

HAPLOTYPES = [1,2] if (config['parameters']['genome_personalization_module']['variant_calling']['make_customRef_diploid'] and (creating_custom_genome or continuing_after_genome_personalization)) else [1] # haplotype number determines number of parallelized runs; if not using a custom PG2-created genome, then num haplotypes = 1

matched_tumor_normal_data=config['user_defined_workflow']['genome_personalization_module']['data_is_matched_tumor_normal']

# inialize the 'experiment/tumor' and 'control/normal' groups
exp_group='tumor' if matched_tumor_normal_data else 'experiment'
ctrl_group='normal' if matched_tumor_normal_data else 'control'

# The transcriptome, genome annotation, and/or gene fusion tracks are merged at the end to comprise the proteome

ALL_TRACKS=[]
BEDFILE_TRACKS=[]

if RNA_seq_module_enabled:
    RNAseq_file_format = config['input_files']['RNA-seq_module']['input_file_format']
    BAM_SAMPLES=[]
    FASTQ_SAMPLES=[]
    EXPERIMENT_SAMPLES=[]
    CONTROL_SAMPLES=[]
    SAMPLE_DICT=dict()
    SAMPLE_REPLICATE_DICT=dict()
    REPLICATE_FILE_DICT=dict()
    SAMPLE_TO_FORMAT=dict()
    if 'bam' in RNAseq_file_format: 
        try: 
            BAM_SAMPLES = list(config['input_files']['RNA-seq_module']['bam_inputs'].keys())
        except KeyError: 
            print("bam specified as an RNA input format, but no bam samples were provided!")
            raise
        for s in BAM_SAMPLES:
            SAMPLE_TO_FORMAT[s] = 'bam'
            
            SAMPLE_REPLICATE_DICT[s] = 'merged'

            if config['input_files']['RNA-seq_module']['bam_inputs'][s]['experiment_vs_control']=='experiment': EXPERIMENT_SAMPLES.append(s)
            elif config['input_files']['RNA-seq_module']['bam_inputs'][s]['experiment_vs_control']=='control': CONTROL_SAMPLES.append(s)
        SAMPLE_DICT[('bam',exp_group)]=EXPERIMENT_SAMPLES
        SAMPLE_DICT[('bam',ctrl_group)]=CONTROL_SAMPLES

    if 'fastq' in RNAseq_file_format:
        try: 
            FASTQ_SAMPLES = list(config['input_files']['RNA-seq_module']['fastq_inputs'].keys())
        except KeyError: 
            print("fastq specified as an RNA input format, but no fastq samples were provided!")
            raise
        for s in FASTQ_SAMPLES:
            SAMPLE_TO_FORMAT[s] = 'fastq'
            
            replicates = config['input_files']['RNA-seq_module']['fastq_inputs'][s]['read_groups'].keys()
            rep_dict_update = {(s,r): (config['input_files']['RNA-seq_module']['fastq_inputs'][s]['read_groups'][r]['R1_fq.gz'],config['input_files']['RNA-seq_module']['fastq_inputs'][s]['read_groups'][r]['R2_fq.gz']) for r in replicates}
            REPLICATE_FILE_DICT.update(rep_dict_update)

            if config['input_files']['RNA-seq_module']['fastq_inputs'][s]['experiment_vs_control']=='experiment': 
                EXPERIMENT_SAMPLES.append(s)
                SAMPLE_REPLICATE_DICT[(s,exp_group)] = replicates
            elif config['input_files']['RNA-seq_module']['fastq_inputs'][s]['experiment_vs_control']=='control': 
                CONTROL_SAMPLES.append(s)
                SAMPLE_REPLICATE_DICT[(s,ctrl_group)] = replicates
        SAMPLE_DICT[('fastq',exp_group)]=list(set(EXPERIMENT_SAMPLES)-set(BAM_SAMPLES))
	SAMPLE_DICT[('fastq',ctrl_group)]=list(set(CONTROL_SAMPLES)-set(BAM_SAMPLES))
    
    STUDY_GROUPS = [exp_group, ctrl_group] if len(CONTROL_SAMPLES)>0 else [exp_group]
    SAMPLE_DICT[exp_group] = EXPERIMENT_SAMPLES
    SAMPLE_DICT[ctrl_group] = CONTROL_SAMPLES
    
    if config['user_defined_workflow']['RNA-seq_module']['transcriptome_track']['assemble_transcriptome_with_StringTie']:
        ALL_TRACKS.append('transcriptome')
        BEDFILE_TRACKS.append('transcriptome')
    if config['user_defined_workflow']['RNA-seq_module']['gene_fusion_track']['assemble_gene_fusions_with_Arriba']:
        ALL_TRACKS.append('gene_fusions')

if config['user_defined_workflow']['genome_annotation_track']['track_enabled']:
    separate_tumor_normal_genomes = config['user_defined_workflow']['genome_personalization_module']['create_separate_tumor_normal_genomes']
    STUDY_GROUPS= [exp_group,ctrl_group] if matched_tumor_normal_data and separate_tumor_normal_genomes else [exp_group]
    ALL_TRACKS.append('genome')
    BEDFILE_TRACKS.append('genome')
    

### End Workflow Control ###


### Input/Output path resolution utils ###

if creating_custom_genome or continuing_after_genome_personalization:
    if config['user_defined_workflow']['genome_personalization_module']['create_separate_tumor_normal_genomes']:
        PG2_GENOME_FASTA = "out/custom_ref/{}.{{study_group}}.H{{htype}}.fa".format(COHORT)
        PG2_GENOME_GTF = "out/custom_ref/{}.{{study_group}}.H{{htype}}.gtf".format(COHORT)
        PG2_STAR_INDEX = "out/custom_ref/{}.{{study_group}}.h-{{htype}}.STARindex/SA".format(COHORT)
        PG2_GENOME_CHAIN = "out/custom_ref/{}.{{study_group}}.H{{htype}}.chain".format(COHORT)
    else:
        PG2_GENOME_FASTA = "out/custom_ref/{}.{}.H{{htype}}.fa".format(COHORT,exp_group)
        PG2_GENOME_GTF = "out/custom_ref/{}.{}.H{{htype}}.gtf".format(COHORT,exp_group)
        PG2_STAR_INDEX = "out/custom_ref/{}.{}.h-{{htype}}.STARindex/SA".format(COHORT,exp_group)
        PG2_GENOME_CHAIN = "out/custom_ref/{}.{}.H{{htype}}.chain".format(COHORT,exp_group)

else:
    PG2_GENOME_FASTA = STOCK_GENOME_FASTA
    PG2_GENOME_GTF = STOCK_GENOME_GTF
    prebuilt_STAR_index_dir = config['stock_references']['genome']['optional_aligner_indices']['STAR_index_dir']
    PG2_STAR_INDEX = os.path.join(prebuilt_STAR_index_dir,'SA') if prebuilt_STAR_index_dir else "out/custom_ref/{}.{{study_group}}.h-{{htype}}.STARindex/SA".format(os.path.basename(PG2_GENOME_FASTA).strip('.fa'))

CHROMOSOMES=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
 
snakemake.utils.makedirs('out/benchmarks')
snakemake.utils.makedirs('out/logs/chr-wise')
snakemake.utils.makedirs('out/logs/RNAseq')
snakemake.utils.makedirs('out/logs/proteomics')

### End path utils ###


### SNAKEMAKE RULES ###

# The 'inputs' of rule 'all' are the default final 'outputs' of PG2; i.e. Snakemake terminates when these files are present
"""
*proteome generation endpoints*
out/{study_group}/combined.proteome.unique.fasta : proteome database FASTA comprising the non-redundant superset of the protein sequences generated by all haplotypes & tracks; used as the input DB during the MaxQuant step
out/{study_group}/combined.proteome.bed : IGV-compatible BED file representation of the aforementioned proteome database

*Proteomics (MaxQuant) endpoint*
out/{study_group}/MaxQuant/combined/txt/summary.txt : MaxQuant (MQ) summary stats file generated at the end of a MQ run; the other MQ output files are also contained in this folder

*Novel peptide analysis endpoints*
out/{study_group}/novel_analysis/{mutation_type}/combined.{mutation_type}.map : tabular map of all mutations of a given type, for which there is peptide-level evidence 
"""
rule all:
    input: expand("out/{study_group}/combined.proteome.unique.fasta",study_group=STUDY_GROUPS), \
           expand("out/{study_group}/combined.proteome.bed",study_group=STUDY_GROUPS), \
           expand("out/{study_group}/novel_analysis/{mutation_type}/{chr}.{mutation_type}.analysis",study_group=STUDY_GROUPS,mutation_type=['missense','insertions','deletions','frameshifts'],chr=CHROMOSOMES), \
           #expand("out/{study_group}/MaxQuant/combined/txt/summary.txt",study_group=STUDY_GROUPS), \
           expand("out/{study_group}/novel_analysis/{mutation_type}/combined.{mutation_type}.map",study_group=STUDY_GROUPS,mutation_type=['missense','insertions','deletions','frameshifts'])
"""
           expand("out/{study_group}/novel_analysis/frameshifts/combined.frameshifts.map",study_group=STUDY_GROUPS), \
           expand("out/{study_group}/novel_analysis/missense/combined.missense.map",study_group=STUDY_GROUPS), \
           expand("out/{study_group}/novel_analysis/insertions/combined.insertions.map",study_group=STUDY_GROUPS), \
           expand("out/{study_group}/novel_analysis/deletions/combined.deletions.map",study_group=STUDY_GROUPS)
"""
# Subworkflows are invoked on rule inputs, and are executed first
subworkflow create_custom_genome:
    snakefile: "genome_personalization_module.py"
    configfile: workflow.overwrite_configfile
    workdir: WD

include: 'novel_analysis.py'

### RNA-seq Mapping & Transcriptome Assembly Workflow ###

rna_seq_read_length = config['input_files']['RNA-seq_module']['read_length']
reads_are_paired_end = config['input_files']['RNA-seq_module']['reads_are_paired_end']

SJ_OVERHANG = int(rna_seq_read_length/2)-1 if reads_are_paired_end else rna_seq_read_length-1
MAX_SEEDS_PER_WINDOW=config['parameters']['RNA-seq_module']['STAR_alignment']['advanced']['max_seeds_per_window']
MAX_INTRON_LENGTH=config['parameters']['RNA-seq_module']['STAR_alignment']['advanced']['max_intron_length']
MAX_MATES_GAP=config['parameters']['RNA-seq_module']['STAR_alignment']['advanced']['max_mates_gap']

rule RNA_00_STAR_CreateGenomeIndex:
    input: fasta=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA),gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
    output: expand("out/custom_ref/{cohort}.{{study_group}}.h-{{htype}}.STARindex/SA",cohort=(COHORT if creating_custom_genome or continuing_after_genome_personalization else os.path.basename(PG2_GENOME_FASTA).strip('.fa')))
    log: "out/logs/{study_group}/h-{htype}.index.txt"
    conda: "envs/STAR.yaml"
    params: directory=os.path.dirname(PG2_STAR_INDEX), n="16", mem_per_cpu="6", R="'span[hosts=1] rusage[mem=6]'", J="index", o="out/logs/index.out", eo="out/logs/index.err"
    shell: "mkdir -p {params.directory} ; \
            STAR \
            --runThreadN {params.n} \
            --runMode genomeGenerate \
            --genomeDir {params.directory} --sjdbGTFfile {input.gtf} --sjdbOverhang {SJ_OVERHANG} --genomeSuffixLengthMax 1000 \
            --genomeFastaFiles {input.fasta} 2> {log}"

if 'bam' in RNAseq_file_format:
    rule RNA_00_ExtractFastqReadsFromRNAseqBAM:
        input: lambda wildcards: config['input_files']['RNA-seq_module']['bam_inputs'][wildcards.sample]['bam_file']
        output: read_one=temp("out/temp_inputs/{sample}.bam2fq.1.fq.gz"),read_two=temp("out/temp_inputs/{sample}.bam2fq.2.fq.gz")
        conda: "envs/STAR.yaml"
        params: n="16", mem_per_cpu="6", R="'span[hosts=1] rusage[mem=6]'", J="RNAseq_bam2fq", o="out/logs/RNAseq/bam2fq.out", eo="out/logs/RNAseq/bam2fq.err",int_readOne=os.path.join(TMP,"{sample}.RG.bam2fq.1.fq"),int_readTwo=os.path.join(TMP,"{sample}.RG.bam2fq.2.fq")
        shell: "samtools collate -O -@ {params.n} {input} | samtools fastq -@ {params.n} -1 {params.int_readOne} -2 {params.int_readTwo} -; gzip -c {params.int_readOne} > {output.read_one}; gzip -c {params.int_readTwo} > {output.read_two}"
 
if 'fastq' not in RNAseq_file_format:
    for group in STUDY_GROUPS:
        SAMPLE_DICT[('fastq',group)] = []
 
if RNA_seq_module_enabled:
    import uuid
    rule RNA_01_STAR_AlignRNAReadsByRG:
        input: PG2_STAR_INDEX, \
               r1 = lambda wildcards: [config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'][replicate]['R1_fq.gz'] for replicate in SAMPLE_REPLICATE_DICT[(wildcards.sample,wildcards.study_group)]] if ('fastq' in RNAseq_file_format and wildcards.sample in SAMPLE_DICT[('fastq',wildcards.study_group)]) else "out/temp_inputs/{sample}.bam2fq.1.fq.gz", \
               r2 = lambda wildcards: [config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'][replicate]['R2_fq.gz'] for replicate in SAMPLE_REPLICATE_DICT[(wildcards.sample,wildcards.study_group)]] if ('fastq' in RNAseq_file_format and wildcards.sample in SAMPLE_DICT[('fastq',wildcards.study_group)]) else "out/temp_inputs/{sample}.bam2fq.2.fq.gz", \
               gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
        output: "out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.sortedByCoord.out.bam"
        conda: "envs/STAR.yaml"
        params: directory=os.path.dirname(PG2_STAR_INDEX), n="16", R="'span[hosts=1] rusage[mem=12]'", J="STAR_align", o="out/logs/RNAseq/{study_group}.haplotype-{htype}.{sample}.STAR.out", eo="out/logs/RNAseq/{study_group}.haplotype-{htype}.{sample}.STAR.err", mem_per_cpu="12", \
                r1_formatted=lambda wildcards:','.join([REPLICATE_FILE_DICT[(wildcards.sample,rep)][0] for rep in SAMPLE_REPLICATE_DICT[(wildcards.sample,wildcards.study_group)]]) if wildcards.sample in SAMPLE_DICT[('fastq',wildcards.study_group)] else "out/temp_inputs/{sample}.bam2fq.1.fq.gz", r2_formatted=lambda wildcards:','.join([REPLICATE_FILE_DICT[(wildcards.sample,rep)][1] for rep in SAMPLE_REPLICATE_DICT[(wildcards.sample,wildcards.study_group)]]) if wildcards.sample in SAMPLE_DICT[('fastq',wildcards.study_group)] else "out/temp_inputs/{sample}.bam2fq.2.fq.gz", \
                tmp_dir=lambda wildcards: os.path.join(TMP,'{}.{}'.format(wildcards.sample,uuid.uuid4()))
        shell: "STAR \
            --genomeDir {params.directory} \
            --readFilesIn {params.r1_formatted} {params.r2_formatted} \
            --outFileNamePrefix out/{wildcards.study_group}/haplotype-{wildcards.htype}/RNAseq/alignment/{wildcards.sample}. \
            --outSAMattributes NH HI XS \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
            --runThreadN {params.n} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMprimaryFlag AllBestScore \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --seedPerWindowNmax {MAX_SEEDS_PER_WINDOW} \
            --outFilterIntronMotifs None \
            --outReadsUnmapped None \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --chimScoreMin 1 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimOutJunctionFormat 1 \
            --chimOutType WithinBAM SoftClip \
            --alignMatesGapMax {MAX_MATES_GAP} \
            --alignIntronMax {MAX_INTRON_LENGTH} \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outFilterType Normal \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 8 \
            --outFilterMismatchNmax 1 \
            --outSJfilterReads Unique \
            --sjdbOverhang {SJ_OVERHANG} \
            --sjdbGTFfile {input.gtf} "
"""
        shell: "STAR \
            --genomeDir {params.directory} \
            --readFilesIn {params.r1_formatted} {params.r2_formatted} \
            --outFileNamePrefix out/{wildcards.study_group}/haplotype-{wildcards.htype}/RNAseq/alignment/{wildcards.sample}. \
            --outTmpDir {params.tmp_dir} \
            --outSAMattributes NH HI XS \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
            --runThreadN {params.n} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMprimaryFlag AllBestScore \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --seedPerWindowNmax {MAX_SEEDS_PER_WINDOW} \
            --outFilterIntronMotifs None \
            --outReadsUnmapped None \
            --chimSegmentMin 10 \
            --chimJunctionOverhangMin 10 \
            --chimScoreMin 1 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimOutJunctionFormat 1 \
            --chimOutType WithinBAM SoftClip \
            --alignMatesGapMax {MAX_MATES_GAP} \
            --alignIntronMax {MAX_INTRON_LENGTH} \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outFilterType Normal \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 8 \
            --outFilterMismatchNmax 1 \
            --outSJfilterReads Unique \
            --sjdbOverhang {SJ_OVERHANG} \
            --sjdbGTFfile {input.gtf} "
"""
# Filter aligned reads in accordance with best practices
import math
max_allowed_multimaps = config['parameters']['RNA-seq_module']['read_filtering']['maximum_allowed_multimaps']
bamflag_filters = config['parameters']['RNA-seq_module']['read_filtering']['advanced']['filter_out_bamFlags']
MIN_MAPQ = (255 if max_allowed_multimaps==1 else int(-10*math.log(1-(1/max_allowed_multimaps),10)))

rule RNA_02_FilterLowQualityReads:
    input: bam="out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.sortedByCoord.out.bam"
    output: "out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.trimmed.out.bam"
    conda: "envs/STAR.yaml"
    params: n="8", mem_per_cpu="2", R="'span[hosts=1] rusage[mem=2]'", J="filter", o="out/logs/filter.out", eo="out/logs/filter.err", flags=bamflag_filters 
    shell: "samtools view -b -h -@ {params.n}\
                $(echo {params.flags} | sed -r 's/[^ ]+/-F &/g') \
                -q {MIN_MAPQ} \
                {input.bam} > {output}"

rule RNA_03_IndexBAMPerSample:
    input: "out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.trimmed.out.bam"
    output: "out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.trimmed.out.bai"
    conda: "envs/STAR.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="BuildBamIndex", o="out/logs/BuildBamIndex.out", eo="out/logs/BuildBamIndex.err"
    shell: "picard BuildBamIndex INPUT={input}"


# Transcriptome Assembly w/ StringTie #

transcriptome_assembly_mode = config['user_defined_workflow']['RNA-seq_module']['transcriptome_track']['GTF-guided-mapping_or_denovo-assembly']

MIN_COVERAGE = config['parameters']['RNA-seq_module']['StringTie_transcriptome_assembly']['advanced']['min_coverage']
MIN_ISOFORM_FRAC = config['parameters']['RNA-seq_module']['StringTie_transcriptome_assembly']['advanced']['min_isoform_frac']

MIN_TPM_MERGE = config['parameters']['RNA-seq_module']['StringTie_transcriptome_assembly']['advanced']['min_tpm_merge']
MIN_FPKM_MERGE = config['parameters']['RNA-seq_module']['StringTie_transcriptome_assembly']['advanced']['min_fpkm_merge']
keeping_retained_introns = config['parameters']['RNA-seq_module']['StringTie_transcriptome_assembly']['advanced']['keep_retained_introns']
# even when mode == 'denovo', guided will run in order to generate the set of fully covered transcripts
rule RNA_04_trscrpt_AssembleWithStringTie_guided:
    input: bam="out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.trimmed.out.bam", bai="out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.trimmed.out.bai",gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
    output: transcriptome="out/{study_group}/haplotype-{htype}/transcriptome/{sample}.StringTie_guided.gtf",covered_refs="out/{study_group}/haplotype-{htype}/transcriptome/{sample}.StringTie_covRefs.gtf"
    conda: "envs/stringtie.yaml"
    params: n="8", mem_per_cpu="4", R="'span[hosts=1]'", J="StringTie", o="out/logs/{study_group}/logsTie_guided.out", eo="out/logs/{study_group}/logsTie_guided.err"
    shell: "stringtie {input.bam} -p {params.n} -o {output.transcriptome} \
                  -G {input.gtf} -C {output.covered_refs} \
                  -c {MIN_COVERAGE} -m {nuc_ORF} -f {MIN_ISOFORM_FRAC}"

if RNA_seq_module_enabled and transcriptome_assembly_mode == 'denovo':
    rule RNA_05_trscrpt_AssembleWithStringTie_denovo:
        input: bam="out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.trimmed.out.bam", bai="out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.trimmed.out.bai"
        output: transcriptome="out/{study_group}/haplotype-{htype}/transcriptome/{sample}.StringTie_denovo.gtf"
        conda: "envs/stringtie.yaml"
        params: n="8", mem_per_cpu="4", R="'span[hosts=1]'", J="StringTie", o="out/logs/{study_group}/logsTie_denovo.out", eo="out/logs/{study_group}/logsTie_denovo.err"
        shell: "stringtie {input.bam} -p {params.n} -o {output.transcriptome} \
                  -c {MIN_COVERAGE} -m {nuc_ORF} -f {MIN_ISOFORM_FRAC}"

if RNA_seq_module_enabled:
    rule RNA_06_trscrpt_CreateSubsetOfFullyCoveredRefTranscripts:
        input: gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF),covered_refs=lambda wildcards: expand("out/{{study_group}}/haplotype-{{htype}}/transcriptome/{sample}.StringTie_covRefs.gtf",sample=SAMPLE_DICT[wildcards.study_group])
        output: gtf_subset=temp("out/{study_group}/haplotype-{htype}/transcriptome/gtf_subset.covRefsOnly.gtf")
        params: n="1", mem_per_cpu="16", R="'span[hosts=1] rusage[mem=16]'", J="subset_refGTF", o="out/logs/subset_refGTF.out", eo="out/logs/subset_refGTF.err"
        shell: "python3 {PG2_HOME}/scripts/subset_fully_covered_transcripts.py {output.gtf_subset} {input.gtf} {input.covered_refs}"

    retaining_fully_covered_refTranscripts = config['user_defined_workflow']['RNA-seq_module']['transcriptome_track']['retain_all_fully_covered_reference_transcripts']
    if retaining_fully_covered_refTranscripts:
        rule RNA_07_trscrpt_MergeSampleWiseTranscriptomesPlusCoveredRefs:
            input: sample_transcriptome=lambda wildcards: expand("out/{{study_group}}/haplotype-{{htype}}/transcriptome/{sample}.StringTie_{mode}.gtf",sample=SAMPLE_DICT[wildcards.study_group],mode=transcriptome_assembly_mode), gtf_subset="out/{study_group}/haplotype-{htype}/transcriptome/gtf_subset.covRefsOnly.gtf"
            output: "out/{study_group}/haplotype-{htype}/transcriptome/transcriptome.gtf"
            log: "out/logs/{study_group}/h-{htype}.merge.txt"
            conda: "envs/stringtie.yaml"
            params: n="8", mem_per_cpu="4", R="'span[hosts=1]'", J="merge", o="out/logs/merge.out", eo="out/logs/merge.err", keep_introns='-i' if keeping_retained_introns else '' 
            shell: "stringtie --merge -o {output} -p {params.n} \
                    -c {MIN_COVERAGE} -m {nuc_ORF} -T {MIN_TPM_MERGE} -F {MIN_FPKM_MERGE} -f {MIN_ISOFORM_FRAC} {params.keep_introns} \
                    -G {input.gtf_subset} \
                    {input.sample_transcriptome} 2> {log}"
    else:
        rule RNA_07_trscrpt_MergeSampleWiseTranscriptomes:
            input: sample_transcriptome=lambda wildcards: expand("out/{{study_group}}/haplotype-{{htype}}/transcriptome/{sample}.StringTie_{mode}.gtf",sample=SAMPLE_DICT[wildcards.study_group],mode=transcriptome_assembly_mode)
            output: "out/{study_group}/haplotype-{htype}/transcriptome/transcriptome.gtf"
            log: "out/logs/{study_group}/h-{htype}.merge.txt"
            conda: "envs/stringtie.yaml"
            params: n="8", mem_per_cpu="4", R="'span[hosts=1]'", J="merge", o="out/logs/merge.out", eo="out/logs/merge.err", keep_introns='-i' if keeping_retained_introns else ''
            shell: "stringtie --merge -o {output} -p {params.n} \
                    -c {MIN_COVERAGE} -m {nuc_ORF} -T {MIN_TPM_MERGE} -F {MIN_FPKM_MERGE} -f {MIN_ISOFORM_FRAC} {params.keep_introns} \
                    {input.sample_transcriptome} 2> {log}"

### END RNA-seq MODULE ###


#TODO: function to create subsets of the genome GTF.
if 'genome' in ALL_TRACKS:
    rule GTF_CreateGenomeAnnotationTrack:
        input: gtf=os.path.abspath((create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF))
        output: "out/{study_group}/haplotype-{htype}/genome/genome.gtf"
        params: n="8", mem_per_cpu="4", R="'span[hosts=1]'", J="merge", o="out/logs/merge.out", eo="out/logs/merge.err"
        shell: "ln -s {input} {output}"


### Gene fusions ###
snakemake.utils.makedirs('out/logs/fusions')

rule arriba:
    input: bam="out/{study_group}/haplotype-{htype}/RNAseq/alignment/{sample}.Aligned.sortedByCoord.out.bam", fa=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA),gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
    output: fusions="out/{study_group}/haplotype-{htype}/gene_fusions/{sample}.fusions.tsv",discarded="out/{study_group}/haplotype-{htype}/gene_fusions/{sample}.fusions_discarded.tsv"
    params: n="8", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", J="arriba_fusions", o="out/logs/fusions/arriba.out", eo="out/logs/arriba/fusions.err", \
            chroms=CHROMOSOMES
    conda: "envs/arriba.yaml"
    #shell: "arriba -x {input.bam} -o {output.fusions} -O {output.discarded} -a {input.fa} -g {input.gtf} -f blacklist"
    #shell: "arriba -x {input.bam} -o {output.fusions} -O {output.discarded} -a {input.fa} -g {input.gtf} -f blacklist -T -P"
    shell: "/home/kwokn/arriba_v1.2.0/arriba -x {input.bam} -o {output.fusions} -O {output.discarded} -a {input.fa} -g {input.gtf} -f blacklist -T -P"

rule CompileFusionTranscripts:
    input: lambda wildcards:expand("out/{{study_group}}/haplotype-{{htype}}/gene_fusions/{sample}.fusions.tsv",sample=SAMPLE_DICT[wildcards.study_group])
    output: "out/{study_group}/haplotype-{htype}/gene_fusions/transcripts.fasta"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", J="compile_fusions", o="out/logs/fusions/compile_fusion_transcripts.out", eo="out/logs/fusions/compile_fusion_transcripts.err"
    shell: "python3 {PG2_HOME}/scripts/compile_fusion_cDNA.py {output} {input}"

rule FormatFusionOrfDB:
    input: orfs="out/{study_group}/haplotype-{htype}/gene_fusions/transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: fasta="out/{study_group}/haplotype-{htype}/gene_fusions/proteome.fasta"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", J="reformat_fusion_ORFs", o="out/logs/fusions/reformat_fusion_ORFs.out", eo="out/logs/fusions/reformat_fusion_ORFs.err"
    run: 
        import textwrap
        with open(input.orfs) as oldDB, open(output.fasta,'w') as newDB:
            for line in oldDB:
                if '>' in line: newDB.write('>{}'.format(line.split(' ')[-1]))
                else: newDB.write('{}\n'.format('\n'.join(textwrap.wrap(line.replace('*','', 80)))))


### Proteome Generation Workflow ###

# Using assembled/annotated exon boundaries, read out the corresponding nucleotide sequences
rule main_01_ExtractCdnaSequences:
    input: gtf="out/{study_group}/haplotype-{htype}/{track}/{track}.gtf", ref_fasta=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA)
    output: fasta = "out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta", gtf="out/{study_group}/haplotype-{htype}/{track}/transcripts.gtf"
    conda: "envs/transdecoder.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="gtf_file_to_cDNA_seqs", o="out/logs/gtf_file_to_cDNA_seqs.out", eo="out/logs/gtf_file_to_cDNA_seqs.err"
    shell: "gffread {input.gtf} -T -o {output.gtf} \
        --force-exons \
        -M -Q; \
        gffread -w {output.fasta} -g {input.ref_fasta} {output.gtf}"

rule main_01a_GTFtoAlignmentGFF3:
    input: "out/{study_group}/haplotype-{htype}/{track}/transcripts.gtf"
    output: "out/{study_group}/haplotype-{htype}/{track}/transcripts.gff3"
    log: "out/logs/{study_group}/h-{htype}.{track}.gtf_to_alignment_gff3.txt"
    conda: "envs/transdecoder.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="gtf_to_alignment_gff3", o="out/logs/gtf_to_alignment_gff3.out", eo="out/logs/gtf_to_alignment_gff3.err"
    shell: "perl {PG2_HOME}/utils/transdecoder/util/gtf_to_alignment_gff3.pl {input} > {output} 2> {log}"

rule main_02_ORF_CalculateCandidateORFs:
    input: "out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta"
    output: "out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir/longest_orfs.pep",checkpoint_dir=directory("out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir.__checkpoints_longorfs/")
    conda: "envs/transdecoder.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="LongOrfs", o="out/logs/LongOrfs.out", eo="out/logs/LongOrfs.err"
    shell: "rm -r {output.checkpoint_dir}; cd out/{wildcards.study_group}/haplotype-{wildcards.htype}/{wildcards.track}; \
        TransDecoder.LongOrfs \
        -t transcripts.fasta \
        -m {ORF}"

PGM_DBNAME = os.path.join(os.path.dirname(STOCK_PROTEOME_FASTA),config['stock_references']['proteome']['fasta'])
rule main_02a_ORF_MakeBlastDB:
    input: fasta=STOCK_PROTEOME_FASTA
    output: [PGM_DBNAME+'.pin', PGM_DBNAME+'.phr', PGM_DBNAME+'.psq']
    log: "out/logs/makeblastdb.txt"
    conda: "envs/blast.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="makeblastdb", o="out/logs/makeblastdb.out", eo="out/logs/makeblastdb.err"
    shell: "makeblastdb \
        -in {input.fasta} \
        -dbtype prot 2> {log} \
        -out {PGM_DBNAME}"

rule main_02b_ORF_BLASTpForHomologyScore:
    input: [PGM_DBNAME+'.pin', PGM_DBNAME+'.phr', PGM_DBNAME+'.psq'], pep="out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: "out/{study_group}/haplotype-{htype}/{track}/blastp.outfmt6"
    log: "out/logs/{study_group}/h-{htype}.{track}.blastp.txt"
    conda: "envs/blast.yaml"
    params: n="18", mem_per_cpu="4", R="'span[ptile=18] rusage[mem=4]'", J="blastp", o="out/logs/blastp.out", eo="out/logs/blastp.err"
    shell: "blastp \
        -num_threads {params.n} \
        -query {input.pep}  \
        -db {PGM_DBNAME}  \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue 1e-5 \
        > {output} 2> {log}"
retaining_single_best_only = config['parameters']['protein_prediction_main']['advanced']['single_best_ORF_only']
rule main_03_ORF_PredictCodingRegions:
    input: orfs="out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        fasta="out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta",
        blastp="out/{study_group}/haplotype-{htype}/{track}/blastp.outfmt6"
    output: "out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder.pep",checkpoint_dir=directory("out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir.__checkpoints/"),gff3="out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder.gff3"
    conda: "envs/transdecoder.yaml"
    params: n="1", mem_per_cpu="18", R="'rusage[mem=18]'", J="Predict", o="out/logs/Predict.out", eo="out/logs/Predict.err",single_best_only='--single_best_only' if retaining_single_best_only else ''
    shell: "rm -r {output.checkpoint_dir};cd out/{wildcards.study_group}/haplotype-{wildcards.htype}/{wildcards.track}; {PG2_HOME}/utils/transdecoder/TransDecoder.Predict.PG2 \
        -t transcripts.fasta \
        --retain_long_orfs_length {nuc_ORF} \
        -v \
        --retain_blastp_hits blastp.outfmt6 \
        {params.single_best_only} \
        --max_overlap_pct 100"

rule main_04_GenerateCDSinGenomeCoords:
    input: gff3="out/{study_group}/haplotype-{htype}/{track}/transcripts.gff3",
        fasta_td="out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta",
        gff3_td="out/{study_group}/haplotype-{htype}/{track}/transcripts.fasta.transdecoder.gff3"
    output: "out/{study_group}/haplotype-{htype}/{track}/transcripts.genome.gff3"
    log: "out/logs/{study_group}/h-{htype}.{track}.cdna_alignment_orf_to_genome_orf.txt"
    conda: "envs/transdecoder.yaml"
    params: n="1", mem_per_cpu="16", R="'rusage[mem=16]'", J="cdna_alignment_orf_to_genome_orf", o="out/logs/cdna_alignment_orf_to_genome_orf.out", eo="out/logs/cdna_alignment_orf_to_genome_orf.err"
    shell: "perl {PG2_HOME}/utils/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl {input.gff3_td} {input.gff3} {input.fasta_td} > {output} 2> {log}"

rule main_05_ReadOutProteomeFASTA:
    input: gff3 = "out/{study_group}/haplotype-{htype}/{track}/transcripts.genome.gff3", ref_fasta=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA)
    output: "out/{study_group}/haplotype-{htype}/{track}/proteome.fasta"
    log: "out/logs/{study_group}/h-{htype}.{track}.gff3_file_to_proteins.txt"
    conda: "envs/transdecoder.yaml"
    params: n="2", mem_per_cpu="8", R="'rusage[mem=8]'", J="gff3_file_to_proteins", o="out/logs/gff3_file_to_proteins.out", eo="out/logs/gff3_file_to_proteins.err"
    shell: "cat {input.gff3} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_proteins.pl --gff3 /dev/stdin --fasta {input.ref_fasta} | egrep -o '^[^*]+' > {output} 2> {log}"

"""
rule remove_duplicate_proteome_entries:
    input: "out/{study_group}/haplotype-{htype}/{track}/proteome.fasta"
    output: "out/{study_group}/haplotype-{htype}/{track}/proteome.unique.fasta"
    log: "out/logs/{study_group}/h-{htype}.{track}.reorderFASTA.txt"
    #conda: "envs/myenv.yaml"
    conda: "envs/R.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="reorderFASTA", o="out/logs/reorderFASTA.out", eo="out/logs/reorderFASTA.err", wd=WD
    script: "{PG2_HOME}/scripts/reorderFASTA.R"
"""

rule main_06_MergeAllProteomeTracksAndRemoveDups:
    input: expand("out/{{study_group}}/haplotype-{htype}/{track}/proteome.fasta", htype=HAPLOTYPES, track=ALL_TRACKS)
    output: "out/{study_group}/combined.proteome.unique.fasta"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="combine_fastas", o="out/logs/combine_fastas.out", eo="out/logs/combine_fastas.err", wd=WD
    shell: "python3 {PG2_HOME}/scripts/reorderFASTA_select_BLAST+ENST.py {output} {input}"

"""
rule combine_assembly_tracks:
    input: expand("out/{{study_group}}/haplotype-{htype}/RNAseq/proteome.fasta", htype=HAPLOTYPES)
    output: "out/{study_group}/combined.assembly.proteome.unique.fasta"
    #conda: "envs/myenv.yaml"
    conda: "envs/R.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="combine_fastas", o="out/logs/combine_fastas.out", eo="out/logs/combine_fastas.err", wd=WD
    script:"{PG2_HOME}/scripts/reorderFASTA.R"
"""

### BEDfile Generation (for IGV) Workflow ###

rule gff3_file_to_bed:
    input: "out/{study_group}/haplotype-{htype}/{track}/transcripts.genome.gff3"
    output: "out/{study_group}/haplotype-{htype}/{track}/proteome_preLiftBack.bed"
    log: "out/logs/{study_group}/h-{htype}.{track}.gff3_file_to_bed.txt"
    conda: "envs/transdecoder.yaml"
    params: n="1", mem_per_cpu="8", R="'rusage[mem=8]'", J="gff3_file_to_bed", o="out/logs/gff3_file_to_bed.out", eo="out/logs/gff3_file_to_bed.err"
    shell: "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_bed.pl /dev/stdin | tail -n +2 > {output} 2> {log}"


if creating_custom_genome or continuing_after_genome_personalization:
    CHAINSWAP=os.path.join(PG2_HOME, config['non-conda_packages']['chainSwap'])
    rule create_reverse_chains:
        input: (create_custom_genome(PG2_GENOME_CHAIN) if creating_custom_genome else PG2_GENOME_CHAIN)
        output: "out/custom_ref/"+COHORT+".{study_group}_H{htype}.chain.reverse"
        params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="reverse_chains", o="out/logs/{study_group}/h-{htype}.reverse_chains.out", eo="out/logs/{study_group}/h-{htype}.reverse_chains.err"
        shell: "{CHAINSWAP} {input} {output}"
    rule liftOver_bed_coords:
        input: bed="out/{study_group}/haplotype-{htype}/{track}/proteome_preLiftBack.bed", chain="out/custom_ref/"+COHORT+".{study_group}_H{htype}.chain.reverse"
        output: "out/{study_group}/haplotype-{htype}/{track}/proteome.bed"
        #conda: "envs/liftover.yaml"
        params: n="1", mem_per_cpu="8", R="'rusage[mem=8]'", J="liftOver_bed", o="out/logs/{study_group}/h-{htype}.{track}.liftOver_bed.out", eo="out/logs/{study_group}/h-{htype}.{track}.liftOver_bed.err", tmp_bed="out/{study_group}/haplotype-{htype}/{track}/proteome_temp.bed"
        shell: "cat {input.bed} | cut -c 3- > {params.tmp_bed}; {PG2_HOME}/utils/liftOver {params.tmp_bed} {input.chain} {output} {output}.unmapped; rm {params.tmp_bed}"
        #shell: "cat {input.bed} | cut -c 3- > {params.tmp_bed}; liftOver {params.tmp_bed} {input.chain} {output} {output}.unmapped; rm {params.tmp_bed}"
else:
    rule rename_bed:
        input: "out/{study_group}/haplotype-{htype}/{track}/proteome_preLiftBack.bed"
        output: "out/{study_group}/haplotype-{htype}/{track}/proteome.bed"
        params: n="1", mem_per_cpu="8", R="'rusage[mem=8]'", J="rename_bed", o="out/logs/{study_group}/h-{htype}.{track}.rename_bed.out", eo="out/logs/{study_group}/h-{htype}.{track}.rename_bed.err"
        shell: "mv {input} {output}"

rule merge_lifted_bedFiles:
    input: expand("out/{{study_group}}/haplotype-{htype}/{track}/proteome.bed",htype=HAPLOTYPES,track=BEDFILE_TRACKS)
    output: "out/{study_group}/combined.proteome.bed"
    #conda: "envs/myenv.yaml"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="merge_proteome_bed", o="out/logs/merge_proteome_bed.out", eo="out/logs/merge_proteome_bed.err"
    shell: "cat {input} | sort -k1,1 -k2,2n > {output}"


### MaxQuant Workflow ###

EXPERIMENT_RAW_DIR = config['input_files']['proteomics_module']['experiment_LCMS_file_directory']
assert EXPERIMENT_RAW_DIR is not None, "missing LCMS_file_directory!"
CONTROL_RAW_DIR = config['input_files']['proteomics_module']['optional_control_LCMS_file_directory'] or ""
PAR = config['input_files']['proteomics_module']['custom_params_xml'] or PG2_HOME + "/MaxQuant/mqpar_template.xml"
MQ = PG2_HOME + "/MaxQuant/bin/MaxQuantCmd.exe"

RAW_FILE_DICT=dict()
RAW_DIR_DICT=dict()
RAW_FILES=[os.path.join(EXPERIMENT_RAW_DIR,f) for f in os.listdir(EXPERIMENT_RAW_DIR) if f.endswith(".raw")]
E_RAW_FILES=[os.path.join(EXPERIMENT_RAW_DIR,f) for f in os.listdir(EXPERIMENT_RAW_DIR) if f.endswith(".raw")]
#RAW_FILE_DICT['experiment']=E_RAW_FILES
RAW_FILE_DICT[exp_group]=E_RAW_FILES
RAW_DIR_DICT[exp_group]=EXPERIMENT_RAW_DIR
if CONTROL_RAW_DIR: 
    C_RAW_FILES=[os.path.join(CONTROL_RAW_DIR,f) for f in os.listdir(CONTROL_RAW_DIR) if f.endswith(".raw")]
    #RAW_FILE_DICT['control']=C_RAW_FILES
    RAW_FILE_DICT[ctrl_group]=C_RAW_FILES
    RAW_DIR_DICT[ctrl_group]=CONTROL_RAW_DIR

MQ_THREADS=str(len(RAW_FILES)) if len(RAW_FILES) >= 16 else '16'

rule copyRawFiles:
    input: raw= lambda wildcards: os.path.join(RAW_DIR_DICT[wildcards.study_group],wildcards.raw_file)
    output: temp("out/{study_group}/MaxQuant/rawfiles/{raw_file}")
    params: n="1", mem_per_cpu="10", R="'span[hosts=1] rusage[mem=10]'", J="copy_raw", o="out/logs/copy_raw.out", eo="out/logs/copy_raw.err", file_dir=directory("out/{study_group}/MaxQuant/rawfiles")
    shell: "cp {input.raw} {params.file_dir}"
    #shell: "ln -s {input.raw} {params.file_dir}"

rule mqpar_conversion:
    input: fasta="out/{study_group}/combined.proteome.unique.fasta",par=PAR
    output: "out/{study_group}/MaxQuant/analysis_ready.mqpar.xml"
    params: n="1", mem_per_cpu="10", R="'span[hosts=1] rusage[mem=10]'", J="mqpar_conversion", o="out/logs/mqpar_conversion.out", eo="out/logs/mqpar_conversion.err"
    run:
        import os
        with open(input.par) as oldMQPar, open(output[0],"w") as newMQPar:
            param_group_line=False
            param_group_list=[]
            RAW_FILES = RAW_FILE_DICT[wildcards.study_group]
            for line in oldMQPar:
                if '<fastaFilePath>' not in line and '<tempFolder>' not in line and '<fixedCombinedFolder>' not in line and '<numThreads>' not in line and '<string>temp</string>' not in line and '<fixedSearchFolder></fixedSearchFolder>' not in line:
                    newMQPar.write(line)
                if '<FastaFileInfo>' in line:
                    newMQPar.write("<fastaFilePath>" + os.getcwd() + "/"+ input.fasta + "</fastaFilePath>\n")
                if '<maxQuantVersion>' in line:
                    newMQPar.write("<tempFolder>" +  TMP + "</tempFolder>\n")
                if '</fastaFilesFirstSearch>' in line:
                    newMQPar.write("<fixedSearchFolder>" +  os.getcwd() + "/out/{}/MaxQuant/search".format(wildcards.study_group) + "</fixedSearchFolder>\n")
                if '<emailFromAddress>' in line:
                    newMQPar.write("<fixedCombinedFolder>"  + os.getcwd() + "/out/{}/MaxQuant".format(wildcards.study_group) + "</fixedCombinedFolder>\n")
                if '<pluginFolder></pluginFolder>' in line:
                    newMQPar.write("<numThreads>"+ MQ_THREADS +"</numThreads>\n")
                if '<filePaths>' in line:
                    for k in range(len(RAW_FILES)):
                        #newMQPar.write("<string>" + RAW_FILES[k] + "</string>\n")
                        newMQPar.write("<string>{}/rawfiles/{}</string>\n".format(os.path.dirname(os.path.abspath(output[0])),os.path.basename(RAW_FILES[k])))
                if '<experiments>' in line:
                    for k in range(len(RAW_FILES)):
                        newMQPar.write("<string></string>\n")
                if '<fractions>' in line:
                    for k in range(len(RAW_FILES)):
                        newMQPar.write("<short>32767</short>\n")
                if '<ptms>' in line:
                    for k in range(len(RAW_FILES)):
                        newMQPar.write("<boolean>False</boolean>\n")
                if '<paramGroupIndices>' in line:
                    param_groups = config['input_files']['proteomics_module']['paramGroups'].keys()
                    substring_group_dict = {config['input_files']['proteomics_module']['paramGroups'][group]['uniquely_identifying_file_substring']: group for group in param_groups}
                    for k in range(len(RAW_FILES)):
                        if len(param_groups)==1: newMQPar.write("<int>0</int>\n")
                        else: 
                            group_assignments=0
                            print(substring_group_dict)
                            print(os.path.basename(RAW_FILES[k]))
                            for substr in substring_group_dict:
                                if substr in os.path.basename(RAW_FILES[k]): 
                                    group_assignments+=1
                                    newMQPar.write("<int>{}</int>\n".format(substring_group_dict[substr]))
                            assert group_assignments == 1, "File {} does not map to a unique parameter group (assigned groups == {}); please double check the input_files->proteomics_module->paramGroups-><group_#>->uniquely_identifying_file_substring parameter".format(RAW_FILES[k],group_assignments)
                if '<parameterGroup>' in line:
                    param_group_line=True
                if param_group_line:
                    param_group_list.append(line)
                    if '<enzymes>' in line:
                        newMQPar.write("<string>{}</string>\n".format(config['input_files']['proteomics_module']['paramGroups'][0]['protease']))
                    if '<fixedModifications>' in line:
                        [newMQPar.write("<string>{}</string>\n".format(x)) for x in config['input_files']['proteomics_module']['paramGroups'][0]['fixedMods']] 
                    if '<variableModifications>' in line:
                        [newMQPar.write("<string>{}</string>\n".format(x)) for x in config['input_files']['proteomics_module']['paramGroups'][0]['variableMods']] 
                if '</parameterGroup>' in line:
                    param_group_line=False
                    additional_proteases=[]
                    for i in range(1,len(config['input_files']['proteomics_module']['paramGroups'].keys())):
                        for param_line in param_group_list:
                            newMQPar.write(param_line)
                            if '<enzymes>' in param_line:
                                newMQPar.write("<string>{}</string>\n".format(config['input_files']['proteomics_module']['paramGroups'][i]['protease']))
                            if '<fixedModifications>' in param_line:
                                [newMQPar.write("<string>{}</string>\n".format(x)) for x in config['input_files']['proteomics_module']['paramGroups'][i]['fixedMods']]
                            if '<variableModifications>' in param_line:
                                [newMQPar.write("<string>{}</string>\n".format(x)) for x in config['input_files']['proteomics_module']['paramGroups'][i]['variableMods']]

rule maxQuant:
    input: lambda wildcards: expand("out/{study_group}/MaxQuant/rawfiles/{raw_file}",study_group=wildcards.study_group, raw_file=[os.path.basename(x) for x in RAW_FILE_DICT[wildcards.study_group]]), par = "out/{study_group}/MaxQuant/analysis_ready.mqpar.xml",db="out/{study_group}/combined.proteome.unique.fasta"
    output: "out/{study_group}/MaxQuant/combined/txt/summary.txt","out/{study_group}/MaxQuant/combined/txt/peptides.txt","out/{study_group}/MaxQuant/combined/txt/proteinGroups.txt"
    singularity: "docker://mono:6.8.0.123"
    params: n=lambda wildcards: str(max(16,min(24,len(RAW_FILE_DICT[wildcards.study_group])))), J="MQ", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/proteomics/mq.{study_group}.out", eo="out/logs/proteomics/mq.{study_group}.err"
    #params: n=lambda wildcards: str(max(16,min(24,len(RAW_FILE_DICT[wildcards.study_group])))), J="MQ", mem_per_cpu="12", R="'span[hosts=1] rusage[mem=12]'", o="out/logs/proteomics/mq.{study_group}.out", eo="out/logs/proteomics/mq.{study_group}.err"
    shell: "mono {MQ} {input.par}"


