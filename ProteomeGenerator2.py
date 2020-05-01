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

HAPLOTYPES = [1,2] if (config['parameters']['genome_personalization_module']['variant_calling']['make_customRef_diploid'] and (creating_custom_genome or continuing_after_genome_personalization)) else [1] # haplotype number determines number of parallelized runs

# The transcriptome, genome annotation, and/or gene fusion tracks are merged at the end to comprise the proteome
TRACKS=[]
if RNA_seq_module_enabled:
    RNAseq_file_format = config['input_files']['RNA-seq_module']['input_file_format']
    BAM_SAMPLES=[]
    FASTQ_SAMPLES=[]
    if 'bam' in RNAseq_file_format: BAM_SAMPLES=list(config['input_files']['RNA-seq_module']['bam_inputs'].keys())
    if 'fastq' in RNAseq_file_format: FASTQ_SAMPLES=list(config['input_files']['RNA-seq_module']['fastq_inputs'].keys())
    RNA_SAMPLES=BAM_SAMPLES+FASTQ_SAMPLES
    
    if config['user_defined_workflow']['RNA-seq_module']['transcriptome_track']['assemble_transcriptome_with_StringTie']:
        TRACKS.append('transcriptome')
    if config['user_defined_workflow']['RNA-seq_module']['gene_fusion_track']['assemble_fusion_genes_with_STAR-Fusion']:
        TRACKS.append('fusions')
if config['user_defined_workflow']['genome_annotation_track']['track_enabled']:
    TRACKS.append('genome')

### End Workflow Control ###


### Input/Output path resolution utils ###

if creating_custom_genome or continuing_after_genome_personalization:
    PG2_GENOME_FASTA = "out/custom_ref/{}_H{{htype}}.fa".format(COHORT)
    PG2_GENOME_GTF = "out/custom_ref/{}_H{{htype}}.gtf".format(COHORT)
    PG2_STAR_INDEX = "out/custom_ref/{}.h-{{htype}}.STARindex/SA".format(COHORT)
    PG2_GENOME_CHAIN = "out/custom_ref/{}_H{{htype}}.chain".format(COHORT)
    
else:
    PG2_GENOME_FASTA = STOCK_GENOME_FASTA
    PG2_GENOME_GTF = STOCK_GENOME_GTF
    prebuilt_STAR_index_dir = config['stock_references']['genome']['optional_aligner_indices']['STAR_index_dir']
    PG2_STAR_INDEX = os.path.join(prebuilt_STAR_index_dir,'SA') if prebuilt_STAR_index_dir else "out/custom_ref/{}.h-{{htype}}.STARindex/SA".format(os.path.basename(PG2_GENOME_FASTA).strip('.fa'))
 
snakemake.utils.makedirs('out')
snakemake.utils.makedirs('out/benchmarks')
snakemake.utils.makedirs('out/logs/chr-wise')
snakemake.utils.makedirs('out/logs/RNAseq')

### End path utils ###


### SNAKEMAKE RULES ###

# Snakemake terminates when these files are present
rule all:
    input: "out/combined.proteome.unique.fasta", "out/combined.proteome.bed", "out/MaxQuant/combined/txt/summary.txt"

# Subworkflows are invoked on rule inputs, and are executed first
subworkflow create_custom_genome:
    snakefile: "genome_personalization_module.py"
    configfile: workflow.overwrite_configfile
    workdir: WD

### Transcriptome Assembly Workflow ###

rna_seq_read_length = config['input_files']['RNA-seq_module']['read_length']
data_is_paired = config['input_files']['RNA-seq_module']['data_is_paired']
SJ_OVERHANG = int(rna_seq_read_length/2)-1 if data_is_paired else rna_seq_read_length-1

MAX_SEEDS_PER_WINDOW=config['parameters']['RNA-seq_module']['STAR_alignment']['advanced']['max_seeds_per_window']
MAX_INTRON_LENGTH=config['parameters']['RNA-seq_module']['STAR_alignment']['advanced']['max_intron_length']
MAX_MATES_GAP=config['parameters']['RNA-seq_module']['STAR_alignment']['advanced']['max_mates_gap']

rule RNA_00_STAR_CreateGenomeIndex:
    input: fasta=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA),gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
    output: expand("out/custom_ref/{cohort}.h-{{htype}}.STARindex/SA",cohort=(COHORT if creating_custom_genome or continuing_after_genome_personalization else os.path.basename(PG2_GENOME_FASTA).strip('.fa')))
    benchmark: "out/benchmarks/h-{htype}.index.txt"
    log: "out/logs/h-{htype}.index.txt"
    conda: "envs/myenv.yaml"
    params: directory=os.path.dirname(PG2_STAR_INDEX), n="16", R="'span[hosts=1] rusage[mem=6]'", J="index", o="out/logs/index.out", eo="out/logs/index.err"
    shell: "mkdir -p {params.directory} ; \
            STAR \
            --runThreadN {params.n} \
            --runMode genomeGenerate --genomeChrBinNbits 10 \
            --genomeDir {params.directory} --sjdbGTFfile {input.gtf} --sjdbOverhang {SJ_OVERHANG} --genomeSuffixLengthMax 1000 \
            --genomeFastaFiles {input.fasta} 2> {log}"

if RNA_seq_module_enabled and 'bam' in RNAseq_file_format:
    rule RNA_00_ExtractFastqReadsFromRNAseqBAM:
        input: lambda wildcards: config['input_files']['RNA-seq_module']['bam_inputs'][wildcards.sample]['bam_file']
        output: read_one=temp("out/bam_inputs/{sample}.RG.bam2fq.1.fq.gz"),read_two=temp("out/bam_inputs/{sample}.RG.bam2fq.2.fq.gz")
        conda: "envs/myenv.yaml"
        params: n="16", R="'span[hosts=1] rusage[mem=6]'", J="RNAseq_bam2fq", o="out/logs/RNAseq/bam2fq.out", eo="out/logs/RNAseq/bam2fq.err",int_readOne=os.path.join(TMP,"{sample}.RG.bam2fq.1.fq"),int_readTwo=os.path.join(TMP,"{sample}.RG.bam2fq.2.fq")
        shell: "samtools collate -O -@ {params.n} {input} | samtools fastq -@ {params.n} -1 {params.int_readOne} -2 {params.int_readTwo} -; gzip -c {params.int_readOne} > {output.read_one}; gzip -c {params.int_readTwo} > {output.read_two}"

    rule wgs_00_CleanRawFastqsWithFastp:
        input: read_one="out/bam_inputs/{sample}.RG.bam2fq.1.fq.gz",read_two="out/bam_inputs/{sample}.RG.bam2fq.2.fq.gz"
        output: read_one=temp("out/bam_inputs/{sample}.RG.bam2fq.fastp.1.fq.gz"), read_two=temp("out/bam_inputs/{sample}.RG.bam2fq.fastp.2.fq.gz")
        params: n="8", R="'span[hosts=1] rusage[mem=4]'", o="out/logs/fastp.out", eo="out/logs/fastp.err", J="fastp_qc"
        conda: "envs/bwa_picard_samtools.yaml"
        shell: "fastp -i {input.read_one} -o {output.read_one} -I {input.read_two} -O {output.read_two}"

    import uuid
    rule RNA_01_STAR_AlignRNAReadsByRG_BAM:
        input: PG2_STAR_INDEX, read_one ="out/bam_inputs/{sample}.{readgroup}.bam2fq.1.fq.gz", read_two="out/bam_inputs/{sample}.{readgroup}.bam2fq.2.fq.gz", gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
        output: temp("out/haplotype-{htype}/RNAseq/alignment/bam_inputs/{sample}.{readgroup}.Aligned.sortedByCoord.out.bam")
        benchmark: "out/benchmarks/h-{htype}.{sample}.{readgroup}.STAR.json"
        log: "out/logs/h-{htype}.{sample}.{readgroup}.STAR.txt"
        conda: "envs/myenv.yaml"
        params: directory=os.path.dirname(PG2_STAR_INDEX), n="32", R="'span[hosts=1] rusage[mem=4]'", J="STAR_align", o="out/logs/STAR_{sample}_bam.out", eo="out/logs/STAR_{sample}_bam.err", \
                strand_field=lambda wildcards: 'None' if config['input_files']['RNA-seq_module']['bam_inputs'][wildcards.sample]['data_is_stranded'] else 'intronMotif', \
                tmp_dir=lambda wildcards: os.path.join(TMP,'{}.{}.{}'.format(wildcards.readgroup,wildcards.sample,uuid.uuid4()))
        shell: "STAR \
            --genomeDir {params.directory} \
            --readFilesIn {input.read_one} {input.read_two} \
            --outFileNamePrefix out/haplotype-{wildcards.htype}/RNAseq/alignment/bam_inputs/{wildcards.sample}.{wildcards.readgroup}. \
            --outTmpDir {params.tmp_dir} \
            --outSAMattributes NH HI XS \
            --outSAMattrRGline ID:{wildcards.readgroup} SM:{wildcards.sample} \
            --runThreadN {params.n} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMprimaryFlag AllBestScore \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --seedPerWindowNmax {MAX_SEEDS_PER_WINDOW} \
            --outSAMstrandField {params.strand_field} \
            --outFilterIntronMotifs None \
            --outReadsUnmapped None \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --chimOutJunctionFormat 1 \
            --chimOutType WithinBAM \
            --alignMatesGapMax {MAX_MATES_GAP} \
            --alignIntronMax {MAX_INTRON_LENGTH} \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outFilterType Normal \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 8 \
            --outFilterMismatchNmax 1 \
            --outSJfilterReads Unique \
            --sjdbOverhang {SJ_OVERHANG} \
            --sjdbGTFfile {input.gtf} \
            \2 > {log} "

if RNA_seq_module_enabled and 'fastq' in RNAseq_file_format:
    import uuid
    rule RNA_01_STAR_AlignRNAReadsByRG_FQ:
        input: PG2_STAR_INDEX, read_one =lambda wildcards: ([config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R1_fq.gz'], config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R2_fq.gz']] if config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R2_fq.gz'] else config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'][wildcards.readgroup]['R1_fq.gz']), gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
        output: "out/haplotype-{htype}/RNAseq/alignment/fastq_inputs/{sample}.{readgroup}.Aligned.sortedByCoord.out.bam"
        benchmark: "out/benchmarks/h-{htype}.{sample}.{readgroup}.STAR.json"
        log: "out/logs/h-{htype}.{sample}.{readgroup}.STAR.txt"
        conda: "envs/myenv.yaml"
        params: directory=os.path.dirname(PG2_STAR_INDEX), n="32", R="'span[hosts=1] rusage[mem=4]'", J="STAR_align", o="out/logs/STAR_{sample}_fq.out", eo="out/logs/STAR_{sample}_fq.err", \
                strand_field=lambda wildcards: 'None' if config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['data_is_stranded'] else 'intronMotif', \
                tmp_dir=lambda wildcards: os.path.join(TMP,'{}.{}.{}'.format(wildcards.readgroup,wildcards.sample,uuid.uuid4()))
        shell: "STAR \
            --genomeDir {params.directory} \
            --readFilesIn {input.read_one} \
            --outFileNamePrefix out/haplotype-{wildcards.htype}/RNAseq/alignment/fastq_inputs/{wildcards.sample}.{wildcards.readgroup}. \
            --outTmpDir {params.tmp_dir} \
            --outSAMattributes NH HI XS \
            --outSAMattrRGline ID:{wildcards.readgroup} SM:{wildcards.sample} \
            --runThreadN {params.n} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMprimaryFlag AllBestScore \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --seedPerWindowNmax {MAX_SEEDS_PER_WINDOW} \
            --outSAMstrandField {params.strand_field} \
            --outFilterIntronMotifs None \
            --outReadsUnmapped None \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --chimOutJunctionFormat 1 \
            --chimOutType WithinBAM \
            --alignMatesGapMax {MAX_MATES_GAP} \
            --alignIntronMax {MAX_INTRON_LENGTH} \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outFilterType Normal \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 8 \
            --outFilterMismatchNmax 1 \
            --outSJfilterReads Unique \
            --sjdbOverhang {SJ_OVERHANG} \
            --sjdbGTFfile {input.gtf} \
            \2 > {log} "
            #--outFilterMatchNminOverLread 0.33 \
            #--outFilterScoreMinOverLread 0.33 \
            #\2 > {log} "

# Filter aligned reads in accordance with best practices
import math
max_allowed_multimaps = config['parameters']['RNA-seq_module']['read_filtering']['maximum_allowed_multimaps']
bamflag_filters = config['parameters']['RNA-seq_module']['read_filtering']['advanced']['filter_out_bamFlags']
rule RNA_02_FilterLowQualityReads:
    input: bam="out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.{readgroup}.Aligned.sortedByCoord.out.bam"
    output: "out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.{readgroup}.Aligned.trimmed.out.bam"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="filter", o="out/logs/filter.out", eo="out/logs/filter.err", min_mapq=(255 if max_allowed_multimaps==1 else int(-10*math.log(1-(1/max_allowed_multimaps),10))), flags=bamflag_filters 
    shell: "samtools view -b -h \
                $(echo {params.flags} | sed -r 's/[^ ]+/-F &/g') \
                -q {params.min_mapq} \
                {input.bam} > {output}"
    #shell: "samtools view -b -h -F 4 -F 256 -F 512 -q 30 {input.bam} > {output}"

if RNA_seq_module_enabled and 'bam' in RNAseq_file_format:
    rule RNA_03_MergeRGsPerSample_BAM:
        input: bam=lambda wildcards: "out/haplotype-{htype}/RNAseq/alignment/bam_inputs/{sample}.RG.Aligned.trimmed.out.bam"
        output: bam="out/haplotype-{htype}/RNAseq/alignment/bam_inputs/{sample}.Aligned.trimmed.RG-merged.out.bam"
        log: "out/logs/h-{htype}.{sample}.MergeRGs.txt"
        params: n="1", R="'rusage[mem=4]'", J="mergeRGs", o="out/logs/mergeRGs.out", eo="out/logs/mergeRGs.err"
        run: 
            in_path = os.path.abspath(input.bam)
            out_path = os.path.abspath(output.bam)
            shell("ln -s {} {}".format(in_path,out_path))

if RNA_seq_module_enabled and 'fastq' in RNAseq_file_format:
    rule RNA_03_MergeRGsPerSample_FQ:
        input: lambda wildcards: expand("out/haplotype-{{htype}}/RNAseq/alignment/fastq_inputs/{{sample}}.{readgroup}.Aligned.sortedByCoord.out.bam", readgroup=config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'].keys())
        #input: lambda wildcards: expand("out/haplotype-{{htype}}/RNAseq/alignment/fastq_inputs/{{sample}}.{readgroup}.Aligned.trimmed.out.bam", readgroup=config['input_files']['RNA-seq_module']['fastq_inputs'][wildcards.sample]['read_groups'].keys())
        output: "out/haplotype-{htype}/RNAseq/alignment/fastq_inputs/{sample}.Aligned.trimmed.RG-merged.out.bam"
        log: "out/logs/h-{htype}.{sample}.MergeRGs.txt"
        conda: "envs/myenv.yaml"
        params: n="1", R="'rusage[mem=4]'", J="mergeRGs", o="out/logs/mergeRGs.out", eo="out/logs/mergeRGs.err"
        shell: "samtools merge {output} {input}"

rule RNA_04_IndexBAMPerSample:
    input: "out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.Aligned.trimmed.RG-merged.out.bam"
    output: "out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.Aligned.trimmed.RG-merged.out.bai"
    log: "out/logs/h-{htype}.{intype}_{sample}.BuildBamIndex.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="BuildBamIndex", o="out/logs/BuildBamIndex.out", eo="out/logs/BuildBamIndex.err"
    shell: "picard \
            BuildBamIndex \
            INPUT={input} 2> {log}"

transcriptome_assembly_mode = config['user_defined_workflow']['RNA-seq_module']['transcriptome_track']['GTF-guided-mapping_or_denovo-assembly']

# even when mode == 'denovo', guided will run in order to generate the set of fully covered transcripts
rule RNA_05_trscrpt_AssembleWithStringTie_guided:
    input: bam="out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.Aligned.trimmed.RG-merged.out.bam", bai="out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.Aligned.trimmed.RG-merged.out.bai",gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
    output: transcriptome="out/haplotype-{htype}/transcriptome/{intype}/{sample}.StringTie_guided.gtf",covered_refs="out/haplotype-{htype}/transcriptome/{intype}/{sample}.StringTie_covRefs.gtf"
    conda: "envs/stringtie.yaml"
    params: n="8", R="'span[hosts=1]'", J="StringTie", o="out/logs/StringTie_guided.out", eo="out/logs/StringTie_guided.err"
    shell: "stringtie {input.bam} -p {params.n} -o {output.transcriptome} \
                  -G {input.gtf} -C {output.covered_refs} \
                  -c 2.5 -m {nuc_ORF} -f 0.01"

if RNA_seq_module_enabled and transcriptome_assembly_mode == 'denovo':
    rule RNA_05_trscrpt_AssembleWithStringTie_denovo:
        input: bam="out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.Aligned.trimmed.RG-merged.out.bam", bai="out/haplotype-{htype}/RNAseq/alignment/{intype}/{sample}.Aligned.trimmed.RG-merged.out.bai"
        output: transcriptome="out/haplotype-{htype}/transcriptome/{intype}/{sample}.StringTie_denovo.gtf"
        conda: "envs/stringtie.yaml"
        params: n="8", R="'span[hosts=1]'", J="StringTie", o="out/logs/StringTie_denovo.out", eo="out/logs/StringTie_denovo.err"
        shell: "stringtie {input.bam} -p {params.n} -o {output.transcriptome} \
                  -c 2.5 -m {nuc_ORF} -f 0.01"

if RNA_seq_module_enabled:
    rule RNA_06_trscrpt_CreateSubsetOfFullyCoveredRefTranscripts:
        input: gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF),covered_refs=[x for subl in [expand("out/haplotype-{{htype}}/transcriptome/bam_inputs/{sample}.StringTie_covRefs.gtf",sample=BAM_SAMPLES),expand("out/haplotype-{{htype}}/transcriptome/fastq_inputs/{sample}.StringTie_covRefs.gtf",sample=FASTQ_SAMPLES)] for x in subl]
        output: gtf_subset=temp("out/haplotype-{htype}/transcriptome/gtf_subset.covRefsOnly.gtf")
        params: n="1", R="'span[hosts=1] rusage[mem=16]'", J="subset_refGTF", o="out/logs/subset_refGTF.out", eo="out/logs/subset_refGTF.err"
        shell: "python3 {PG2_HOME}/scripts/subset_fully_covered_transcripts.py {output.gtf_subset} {input.gtf} {input.covered_refs}"

    retaining_fully_covered_refTranscripts = config['user_defined_workflow']['RNA-seq_module']['transcriptome_track']['retain_all_fully_covered_reference_transcripts']
    if retaining_fully_covered_refTranscripts:
        rule RNA_07_trscrpt_MergeSampleWiseTranscriptomesPlusCoveredRefs:
            input: sample_transcriptome=[x for subl in [expand("out/haplotype-{{htype}}/transcriptome/bam_inputs/{sample}.StringTie_{mode}.gtf",sample=BAM_SAMPLES,mode=transcriptome_assembly_mode),expand("out/haplotype-{{htype}}/transcriptome/fastq_inputs/{sample}.StringTie_{mode}.gtf",sample=FASTQ_SAMPLES,mode=transcriptome_assembly_mode)] for x in subl], gtf_subset="out/haplotype-{htype}/transcriptome/gtf_subset.covRefsOnly.gtf"
            output: "out/haplotype-{htype}/transcriptome/transcriptome.gtf"
            log: "out/logs/h-{htype}.merge.txt"
            conda: "envs/stringtie.yaml"
            params: n="8", R="'span[hosts=1]'", J="merge", o="out/logs/merge.out", eo="out/logs/merge.err"
            shell: "stringtie --merge -o {output} -p {params.n} \
                    -c 2.5 -m {nuc_ORF} -T 1 -f 0.01 -i \
                    -G {input.gtf_subset} \
                    {input.sample_transcriptome} 2> {log}"
    else:
        rule RNA_07_trscrpt_MergeSampleWiseTranscriptomes:
            input: sample_transcriptome=[x for subl in [expand("out/haplotype-{{htype}}/transcriptome/bam_inputs/{sample}.StringTie_{mode}.gtf",sample=BAM_SAMPLES,mode=transcriptome_assembly_mode),expand("out/haplotype-{{htype}}/transcriptome/fastq_inputs/{sample}.StringTie_{mode}.gtf",sample=FASTQ_SAMPLES,mode=transcriptome_assembly_mode)] for x in subl]
            output: "out/haplotype-{htype}/transcriptome/transcriptome.gtf"
            log: "out/logs/h-{htype}.merge.txt"
            conda: "envs/stringtie.yaml"
            params: n="8", R="'span[hosts=1]'", J="merge", o="out/logs/merge.out", eo="out/logs/merge.err"
            shell: "stringtie --merge -o {output} -p {params.n} \
                    -c 2.5 -m {nuc_ORF} -T 1 -f 0.01 -i \
                    {input.sample_transcriptome} 2> {log}"

#TODO: function to create subsets of the genome GTF.
if 'genome' in TRACKS:
    rule GTF_CreateGenomeAnnotationTrack:
        input: gtf=os.path.abspath((create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF))
        output: "out/haplotype-{htype}/genome/genome.gtf"
        log: "out/logs/h-{htype}.merge.txt"
        params: n="8", R="'span[hosts=1]'", J="merge", o="out/logs/merge.out", eo="out/logs/merge.err"
        shell: "ln -s {input} {output}"

### Gene fusions ###
rule RNA_00_fusion_BuildCTATGenomelib:
    input: fasta=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA),gtf=(create_custom_genome(PG2_GENOME_GTF) if creating_custom_genome else PG2_GENOME_GTF)
    output: "out/custom_ref/haplotype-{htype}/ctat_genome_lib_build_dir/__chkpts/validate_ctat_genome_lib.ok"
    output: "out/custom_ref/haplotype-{htype}/ctat_genome_lib_build_dir/h-{htype}.ref_genome.fa"
    params: n="16", R="'span[hosts=1] rusage[mem=8]'", J="build_fusion_lib", o="out/logs/build_fusion_lib.out", eo="out/logs/build_fusion_lib.err", \
            fasta=os.path.abspath(PG2_GENOME_FASTA),gtf=os.path.abspath(PG2_GENOME_GTF)
    singularity: "docker://trinityctat/starfusion"
    shell: "cd out/custom_ref/haplotype-{wildcards.htype}; perl {PG2_HOME}/utils/CTAT_genome_utils/prep_genome_lib.pl \
                --genome_fa {params.fasta} --gtf {params.gtf} \
                --fusion_annot_lib {PG2_HOME}/utils/CTAT_genome_utils/fusion_annot_lib.gz \
                --annot_filter_rule {PG2_HOME}/utils/CTAT_genome_utils/AnnotFilterRule.pm \
                --pfam_db current --dfam_db human --human_gencode_filter --CPU 8"


### Proteome Generation Workflow ###

# Read out nucleotide sequences from GTFs
rule main_01_ExtractCdnaSequences:
    input: gtf="out/haplotype-{htype}/{track}/{track}.gtf", ref_fasta=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA)
    output: fasta = "out/haplotype-{htype}/{track}/transcripts.fasta",
        gtf="out/haplotype-{htype}/{track}/transcripts.gtf"
    benchmark: "out/benchmarks/h-{htype}.{track}.gtf_file_to_cDNA_seqs.txt"
    log: "out/logs/h-{htype}.{track}.gtf_file_to_cDNA_seqs.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="gtf_file_to_cDNA_seqs", o="out/logs/gtf_file_to_cDNA_seqs.out", eo="out/logs/gtf_file_to_cDNA_seqs.err"
    shell: "gffread {input.gtf} -T -o {output.gtf} \
        --no-pseudo \
        --force-exons \
        -M -Q; \
        gffread -w {output.fasta} -g {input.ref_fasta} {output.gtf} 2> {log}"

rule main_01a_GTFtoAlignmentGFF3:
    input: "out/haplotype-{htype}/{track}/transcripts.gtf"
    output: "out/haplotype-{htype}/{track}/transcripts.gff3"
    benchmark: "out/benchmarks/h-{htype}.{track}.gtf_to_alignment_gff3.txt"
    log: "out/logs/h-{htype}.{track}.gtf_to_alignment_gff3.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="gtf_to_alignment_gff3", o="out/logs/gtf_to_alignment_gff3.out", eo="out/logs/gtf_to_alignment_gff3.err"
    shell: "perl {PG2_HOME}/utils/transdecoder/util/gtf_to_alignment_gff3.pl {input} > {output} 2> {log}"

rule main_02_ORF_CalculateCandidateORFs:
    input: "out/haplotype-{htype}/{track}/transcripts.fasta"
    output: "out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir/longest_orfs.pep",checkpoint_dir=directory("out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir.__checkpoints_longorfs/")
    benchmark: "out/benchmarks/h-{htype}.{track}.LongOrfs.json"
    log: "../../logs/h-{htype}.{track}.LongOrfs.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="LongOrfs", o="out/logs/LongOrfs.out", eo="out/logs/LongOrfs.err"
    shell: "rm -r {output.checkpoint_dir}; cd out/haplotype-{wildcards.htype}/{wildcards.track}; \
        TransDecoder.LongOrfs \
        -t transcripts.fasta \
        -m {ORF} 2> {log}"

PGM_DBNAME = os.path.join(os.path.dirname(STOCK_PROTEOME_FASTA),config['stock_references']['proteome']['fasta'])
rule main_02a_ORF_MakeBlastDB:
    input: fasta=STOCK_PROTEOME_FASTA
    output: [PGM_DBNAME+'.pin', PGM_DBNAME+'.phr', PGM_DBNAME+'.psq']
    benchmark: "out/benchmarks/makeblastdb.json"
    log: "out/logs/makeblastdb.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="makeblastdb", o="out/logs/makeblastdb.out", eo="out/logs/makeblastdb.err"
    shell: "makeblastdb \
        -in {input.fasta} \
        -dbtype prot 2> {log} \
        -out {PGM_DBNAME}"

rule main_02b_ORF_BLASTpForHomologyScore:
    input: [PGM_DBNAME+'.pin', PGM_DBNAME+'.phr', PGM_DBNAME+'.psq'], pep="out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: "out/haplotype-{htype}/{track}/blastp.outfmt6"
    benchmark: "out/benchmarks/h-{htype}.{track}.blastp.json"
    log: "out/logs/h-{htype}.{track}.blastp.txt"
    conda: "envs/myenv.yaml"
    params: n="18", R="'span[ptile=18] rusage[mem=4]'", J="blastp", o="out/logs/blastp.out", eo="out/logs/blastp.err"
    shell: "blastp \
        -num_threads {params.n} \
        -query {input.pep}  \
        -db {PGM_DBNAME}  \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue 1e-5 \
        > {output} 2> {log}"

rule main_03_ORF_PredictCodingRegions:
    input: orfs="out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        fasta="out/haplotype-{htype}/{track}/transcripts.fasta",
        blastp="out/haplotype-{htype}/{track}/blastp.outfmt6"
    output: "out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder.pep",checkpoint_dir=directory("out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder_dir.__checkpoints/"),gff3="out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder.gff3"
    benchmark: "out/benchmarks/h-{htype}.{track}.Predict.json"
    log: "../../logs/h-{htype}.{track}.Predict.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=18]'", J="Predict", o="out/logs/Predict.out", eo="out/logs/Predict.err"
    shell: "rm -r {output.checkpoint_dir};cd out/haplotype-{wildcards.htype}/{wildcards.track}; {PG2_HOME}/utils/transdecoder/TransDecoder.Predict.IGNORE_OVERLAP \
        -t transcripts.fasta \
        --retain_long_orfs_length {nuc_ORF} \
        -v \
        --retain_blastp_hits blastp.outfmt6 2> {log}"

rule main_04_GenerateCDSinGenomeCoords:
    input: gff3="out/haplotype-{htype}/{track}/transcripts.gff3",
        fasta_td="out/haplotype-{htype}/{track}/transcripts.fasta",
        gff3_td="out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder.gff3"
    output: "out/haplotype-{htype}/{track}/transcripts.genome.gff3"
    benchmark: "out/benchmarks/h-{htype}.{track}.cdna_alignment_orf_to_genome_orf.txt"
    log: "out/logs/h-{htype}.{track}.cdna_alignment_orf_to_genome_orf.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=16]'", J="cdna_alignment_orf_to_genome_orf", o="out/logs/cdna_alignment_orf_to_genome_orf.out", eo="out/logs/cdna_alignment_orf_to_genome_orf.err"
    shell: "perl {PG2_HOME}/utils/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl {input.gff3_td} {input.gff3} {input.fasta_td} > {output} 2> {log}"

##TODO: REMOVE THIS & TURN INTO USER-CONTROLLED VARIABLE
incorporatingSomaticVariants=False

if not incorporatingSomaticVariants:
    rule main_05_ReadOutProteomeFASTA:
        input: gff3 = "out/haplotype-{htype}/{track}/transcripts.genome.gff3", ref_fasta=(create_custom_genome(PG2_GENOME_FASTA) if creating_custom_genome else PG2_GENOME_FASTA)
        output: "out/haplotype-{htype}/{track}/proteome.fasta"
        benchmark: "out/benchmarks/h-{htype}.{track}.gff3_file_to_proteins.txt"
        log: "out/logs/h-{htype}.{track}.gff3_file_to_proteins.txt"
        conda: "envs/myenv.yaml"
        params: n="2", R="'rusage[mem=8]'", J="gff3_file_to_proteins", o="out/logs/gff3_file_to_proteins.out", eo="out/logs/gff3_file_to_proteins.err"
        shell: "cat {input.gff3} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_proteins.pl --gff3 /dev/stdin --fasta {input.ref_fasta} | egrep -o '^[^*]+' > {output} 2> {log}"
else:
    rule extract_VCFexpanded_transcripts:
        input: "out/haplotype-{htype}/{track}/transcripts.fasta.transdecoder.gff3"
        output: main="out/haplotype-{htype}/{track}/main_transcripts.gff3",expanded="out/haplotype-{htype}/{track}/expanded_transcripts.gff3"
        log:"out/logs/h-{htype}.{track}.partition_expanded.txt"
        params: n="1", R="'rusage[mem=4]'", J="partition_expanded", o="out/logs/partition_expanded.out", eo="out/logs/partition_expanded.err"
        shell: "python3 {PG2_HOME}/scripts/partition_transcripts_gff3.py {input} {output.main} {output.expanded}"
    rule prune_genome_gff3:
        input: "out/haplotype-{htype}/{track}/transcripts.genome.gff3"
        output: main="out/haplotype-{htype}/{track}/main_transcripts.genome.gff3",expanded="out/haplotype-{htype}/{track}/expanded_transcripts.genome.gff3"
        log:"out/logs/h-{htype}.{track}.prune_genome_gff3.txt"
        params: n="1", R="'rusage[mem=4]'", J="partition_expanded", o="out/logs/partition_expanded.out", eo="out/logs/partition_expanded.err"
        shell: "python3 {PG2_HOME}/scripts/partition_transcripts_gff3.py {input} {output.main} {output.expanded}"
    rule read_out_main_proteome:
        input: gff3 = "out/haplotype-{htype}/{track}/main_transcripts.genome.gff3", ref_fasta=resolve_custom('fa','{htype}')
        output: "out/haplotype-{htype}/{track}/proteome.main.fasta"
        benchmark: "out/benchmarks/h-{htype}.{track}.gff3_file_to_proteins.txt"
        log: "out/logs/h-{htype}.{track}.gff3_file_to_proteins.txt"
        conda: "envs/myenv.yaml"
        params: n="2", R="'rusage[mem=8]'", J="gff3_file_to_proteins", o="out/logs/gff3_file_to_proteins.out", eo="out/logs/gff3_file_to_proteins.err"
        shell: "cat {input.gff3} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_proteins.pl --gff3 /dev/stdin --fasta {input.ref_fasta} | egrep -o '^[^*]+' > {output} 2> {log}"
	
    rule read_out_expanded_transcripts:
        input: gff3 = "out/haplotype-{htype}/{track}/expanded_transcripts.gff3", ref_fasta='out/haplotype-{htype}/{track}/expanded.fasta'
        output: "out/haplotype-{htype}/{track}/proteome.expanded.fasta"
        benchmark: "out/benchmarks/h-{htype}.{track}.gff3_file_to_proteins.txt"
        log: "out/logs/h-{htype}.{track}.gff3_file_to_proteins.txt"
        conda: "envs/myenv.yaml"
        params: n="2", R="'rusage[mem=8]'", J="gff3_file_to_proteins", o="out/logs/gff3_file_to_proteins.out", eo="out/logs/gff3_file_to_proteins.err"
        shell: "samtools faidx {input.ref_fasta}; cat {input.gff3} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_proteins.pl --gff3 /dev/stdin --fasta {input.ref_fasta} | egrep -o '^[^*]+' > {output} 2> {log}"

    rule combine_main_and_expanded:
        input: main="out/haplotype-{htype}/{track}/proteome.main.fasta",expanded="out/haplotype-{htype}/{track}/proteome.expanded.fasta"
        output: "out/haplotype-{htype}/{track}/proteome.fasta"
        params: n="1", R="'rusage[mem=8]'", J="gff3_file_to_proteins", o="out/logs/gff3_file_to_proteins.out", eo="out/logs/gff3_file_to_proteins.err"
        shell: "cat {input.main} {input.expanded} > {output}"
###

rule remove_duplicate_proteome_entries:
    input: "out/haplotype-{htype}/{track}/proteome.fasta"
    output: "out/haplotype-{htype}/{track}/proteome.unique.fasta"
    benchmark: "out/benchmarks/h-{htype}.{track}.reorderFASTA.txt"
    log: "out/logs/h-{htype}.{track}.reorderFASTA.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="reorderFASTA", o="out/logs/reorderFASTA.out", eo="out/logs/reorderFASTA.err", wd=WD
    script: "{PG2_HOME}/scripts/reorderFASTA.R"

rule main_06_MergeAllProteomeTracksAndRemoveDups:
    input: expand("out/haplotype-{htype}/{track}/proteome.fasta", htype=HAPLOTYPES, track=TRACKS)
    output: "out/combined.proteome.unique.fasta"
    benchmark: "out/benchmarks/combine_FASTAs.txt"
    log: "out/logs/combine_FASTAs.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="combine_fastas", o="out/logs/combine_fastas.out", eo="out/logs/combine_fastas.err", wd=WD
    shell: "python3 {PG2_HOME}/scripts/reorderFASTA_select_BLAST+ENST.py {output} {input}"
    #script:"{PG2_HOME}/scripts/reorderFASTA.R"

rule combine_assembly_tracks:
    input: expand("out/haplotype-{htype}/RNAseq/proteome.fasta", htype=HAPLOTYPES)
    output: "out/combined.assembly.proteome.unique.fasta"
    benchmark: "out/benchmarks/combine_FASTAs.txt"
    log: "out/logs/combine_FASTAs.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="combine_fastas", o="out/logs/combine_fastas.out", eo="out/logs/combine_fastas.err", wd=WD
    script:"{PG2_HOME}/scripts/reorderFASTA.R"

rule re_reorderFASTA:
    input: "out/combined.concat_headers.proteome.unique.fasta"
    output: "out/combined.concat_re-reordered.proteome.unique.fasta"
    benchmark: "out/benchmarks/combine_FASTAs.txt"
    log: "out/logs/combine_FASTAs.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="combine_fastas", o="out/logs/combine_fastas.out", eo="out/logs/combine_fastas.err", wd=WD
    script:"{PG2_HOME}/scripts/reorderFASTA.R"


### BEDfile Generation (for IGV) Workflow ###

rule gff3_file_to_bed:
    input: "out/haplotype-{htype}/{track}/transcripts.genome.gff3"
    output: "out/haplotype-{htype}/{track}/proteome_preLiftBack.bed"
    benchmark: "out/benchmarks/h-{htype}.{track}.gff3_file_to_bed.txt"
    log: "out/logs/h-{htype}.{track}.gff3_file_to_bed.txt"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=8]'", J="gff3_file_to_bed", o="out/logs/gff3_file_to_bed.out", eo="out/logs/gff3_file_to_bed.err"
    shell: "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_bed.pl /dev/stdin | tail -n +2 > {output} 2> {log}"


if creating_custom_genome or continuing_after_genome_personalization:
    CHAINSWAP=os.path.join(PG2_HOME, config['non-conda_packages']['chainSwap'])
    rule create_reverse_chains:
        input: (create_custom_genome(PG2_GENOME_CHAIN) if creating_custom_genome else PG2_GENOME_CHAIN)
        output: "out/custom_ref/"+COHORT+"_H{htype}.chain.reverse"
        params: n="1", R="'rusage[mem=4]'", J="reverse_chains", o="out/logs/h-{htype}.reverse_chains.out", eo="out/logs/h-{htype}.reverse_chains.err"
        shell: "{CHAINSWAP} {input} {output}"
    rule liftOver_bed_coords:
        input: bed="out/haplotype-{htype}/{track}/proteome_preLiftBack.bed", chain="out/custom_ref/"+COHORT+"_H{htype}.chain.reverse"
        output: "out/haplotype-{htype}/{track}/proteome.bed"
        conda: "envs/myenv.yaml"
        params: n="1", R="'rusage[mem=8]'", J="liftOver_bed", o="out/logs/h-{htype}.{track}.liftOver_bed.out", eo="out/logs/h-{htype}.{track}.liftOver_bed.err", tmp_bed="out/haplotype-{htype}/{track}/proteome_temp.bed"
        shell: "cat {input.bed} | cut -c 3- > {params.tmp_bed}; liftOver {params.tmp_bed} {input.chain} {output} {output}.unmapped; rm {params.tmp_bed}"
else:
    rule rename_bed:
        input: "out/haplotype-{htype}/{track}/proteome_preLiftBack.bed"
        output: "out/haplotype-{htype}/{track}/proteome.bed"
        params: n="1", R="'rusage[mem=8]'", J="rename_bed", o="out/logs/h-{htype}.{track}.rename_bed.out", eo="out/logs/h-{htype}.{track}.rename_bed.err"
        shell: "mv {input} {output}"

rule merge_lifted_bedFiles:
    input: expand("out/haplotype-{htype}/{track}/proteome.bed",htype=HAPLOTYPES,track=TRACKS)
    output: "out/combined.proteome.bed"
    conda: "envs/myenv.yaml"
    params: n="1", R="'rusage[mem=4]'", J="merge_proteome_bed", o="out/logs/merge_proteome_bed.out", eo="out/logs/merge_proteome_bed.err"
    shell: "cat {input} | sort -k1,1 -k2,2n > {output}"


### MaxQuant Workflow ###

RAW_DIR = config['input_files']['proteomics_module']['LCMS_file_directory']
assert RAW_DIR is not None, "missing LCMS_file_directory!"
PAR = config['input_files']['proteomics_module']['custom_params_xml'] or PG2_HOME + "/MaxQuant/mqpar_template.xml"


RAW_FILES=[f for f in os.listdir(RAW_DIR) if f.endswith(".raw")]

rule copyRawFiles:
    input: raw=os.path.join(RAW_DIR,'{raw_file}'),fasta='out/combined.proteome.unique.fasta'
    output: temp("out/MaxQuant/{raw_file}")
    params: n="1", R="'span[hosts=1] rusage[mem=10]'", J="copy_raw", o="out/logs/copy_raw.out", eo="out/logs/copy_raw.err"
    shell: "cp {input.raw} {output}"

rule mqpar_conversion:
    input: fasta="out/combined.proteome.unique.fasta"
    output: "out/MaxQuant/analysis_ready.mqpar.xml"
    benchmark: "out/benchmarks/mqpar_conversion.txt"
    log: "out/logs/mqpar_conversion.txt"
    params: n="1", R="'span[hosts=1] rusage[mem=10]'", J="mqpar_conversion", o="out/logs/mqpar_conversion.out", eo="out/logs/mqpar_conversion.err"
    run:
        import os
        with open(PAR) as oldMQPar, open(output[0],"w") as newMQPar:
            for line in oldMQPar:
                if '<fastaFilePath>' not in line and '<tempFolder>' not in line and '<fixedCombinedFolder>' not in line and '<numThreads>' not in line and '<string>temp</string>' not in line and '<fixedSearchFolder></fixedSearchFolder>' not in line:
                    newMQPar.write(line)
                if '<FastaFileInfo>' in line:
                    newMQPar.write("<fastaFilePath>" + os.getcwd() + "/"+ input.fasta + "</fastaFilePath>\n")
                if '<maxQuantVersion>' in line:
                    newMQPar.write("<tempFolder>" +  TMP + "</tempFolder>\n")
                if '</fastaFilesFirstSearch>' in line:
                    newMQPar.write("<fixedSearchFolder>" +  os.getcwd() + "/out/MaxQuant/search" + "</fixedSearchFolder>\n")
                if '<emailFromAddress>' in line:
                    newMQPar.write("<fixedCombinedFolder>"  + os.getcwd() + "/out/MaxQuant" + "</fixedCombinedFolder>\n")
                if '<pluginFolder></pluginFolder>' in line:
                    newMQPar.write("<numThreads>"+ THREADS +"</numThreads>\n")
                if '<filePaths>' in line:
                    for k in range(len(RAW_FILES)):
                        newMQPar.write("<string>" + os.getcwd() + "/out/MaxQuant/" + RAW_FILES[k] + "</string>\n")
                if '<experiments>' in line:
                    for k in range(len(RAW_FILES)-1):
                        newMQPar.write("<string></string>\n")
                if '<fractions>' in line:
                    for k in range(len(RAW_FILES)-1):
                        newMQPar.write("<short>32767</short>\n")
                if '<ptms>' in line:
                    for k in range(len(RAW_FILES)-1):
                        newMQPar.write("<boolean>False</boolean>\n")
                if '<paramGroupIndices>' in line:
                    for k in range(len(RAW_FILES)-1):
                        newMQPar.write("<int>0</int>\n")


MQ = PG2_HOME + "/MaxQuant/bin/MaxQuantCmd.exe"
THREADS=str(len(RAW_FILES)) if len(RAW_FILES) >= 16 else '16'
rule maxQuant:
    input: expand("out/MaxQuant/{raw_file}",raw_file=RAW_FILES), par = "out/MaxQuant/analysis_ready.mqpar.xml"
    output: "out/MaxQuant/combined/txt/summary.txt"
    benchmark: "out/benchmarks/maxQuant.txt"
    log: "out/logs/maxQuant.txt"
    singularity: "docker://mono:5.12.0.226"
    params: n=THREADS, J="MQ", R="'span[hosts=1] rusage[mem=8]'".format(THREADS), o="out/logs/mq.out", eo="out/logs/mq.err"
    shell: "mono {MQ} {input.par}"


