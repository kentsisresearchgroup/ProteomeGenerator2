import os
PG2_HOME = config['directories']['PG2_installation_dir']
WD = config['directories']['working_and_output_dir']
workdir: WD

STOCK_GENOME_FASTA=config['stock_references']['genome']['fasta']
STOCK_GENOME_GTF=config['stock_references']['genome']['gtf']


input_file_format = config['input_files']['genome_personalization_module']['input_file_format']
try: [config['input_files']['genome_personalization_module'][fmt+'_inputs'] for fmt in input_file_format]
except: raise TypeError("ERROR: Specified inputs do not match the specified input file format.")

COHORT = config['input_files']['genome_personalization_module']['cohort_or_organism_name']
TUMOR_SAMPLES=[]
NORMAL_SAMPLES=[]
NONMATCHED_SAMPLES=[]
sample_dict_list = [dict(config['input_files']['genome_personalization_module'][fmt+'_inputs']) for fmt in input_file_format]
for d in sample_dict_list:
    for sample_name in d:
        if d[sample_name]['matched_sample_params']['is_matched_sample'] and d[sample_name]['matched_sample_params']['tumor_or_normal']=='tumor': TUMOR_SAMPLES.append(sample_name)
        else:
            NORMAL_SAMPLES.append(sample_name)
            NONMATCHED_SAMPLES.append(sample_name)


ALL_SAMPLES=TUMOR_SAMPLES + NORMAL_SAMPLES

ANALYSIS_READY_BAMFILES=[]
SAMPLEFILE_SAMPLENAME_DICT=dict()
for sample in TUMOR_SAMPLES:
    filename = "out/tumor/variant_calling/{}.analysis_ready.bam".format(sample)
    ANALYSIS_READY_BAMFILES.append(filename)
    SAMPLEFILE_SAMPLENAME_DICT[filename] = sample
for sample in NORMAL_SAMPLES:
    sample_dict=None
    for d in sample_dict_list:
        if sample in d: 
            sample_dict=d
    filename = "out/normal/variant_calling/{}.analysis_ready.bam".format(sample) if sample_dict[sample]['matched_sample_params']['is_matched_sample'] else "out/experiment/variant_calling/{}.analysis_ready.bam".format(sample)
    ANALYSIS_READY_BAMFILES.append(filename)
    SAMPLEFILE_SAMPLENAME_DICT[filename] = sample

VARIANT_CALLING_MODES=[]
if config['user_defined_workflow']['genome_personalization_module']['variant_calling_submodule']['call_germline_variants_with_GATK4_HaplotypeCaller']==True: VARIANT_CALLING_MODES.append('germline')
if config['user_defined_workflow']['genome_personalization_module']['variant_calling_submodule']['call_somatic_variants_with_GATK4_Mutect2']==True: VARIANT_CALLING_MODES.append('somatic')

running_preprocessing = config['user_defined_workflow']['genome_personalization_module']['variant_calling_submodule']['run_pre-processing_steps']
just_ran_WGS_preprocessing = config['user_defined_workflow']['genome_personalization_module']['variant_calling_submodule']['continuation']['just_ran_PG2_preprocessing']

subworkflow WGS_preprocessing:
    snakefile: "WGS_preprocessing.py"
    configfile: workflow.overwrite_configfile
    workdir: WD

snakemake.utils.makedirs('out/logs/variant_calling/intervals')
snakemake.utils.makedirs('out/logs/variant_calling')
snakemake.utils.makedirs('out/logs/variant_calling/chr-wise')
snakemake.utils.makedirs('out/logs/variant_calling/intervals')
NUM_VARIANT_INTERVALS = config['parameters']['genome_personalization_module']['variant_calling']['advanced']['variant_intervals_scatter']
rule var_00_ScatterVariantCallingIntervals:
    input: config['parameters']['genome_personalization_module']['variant_calling']['resources']['wgs_calling_regions']
    output: "out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list"
    params: n="1", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/variant_calling/intervals/{interval}.out", eo="out/logs/variant_calling/intervals/{interval}.err", \
            J="generate_intervals_{interval}"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx8g' SplitIntervals \
              -R {STOCK_GENOME_FASTA} -L {input} -scatter {NUM_VARIANT_INTERVALS} -O out/{wildcards.study_group}/variant_calling/intervals"

if (not just_ran_WGS_preprocessing) and running_preprocessing:
    rule var_00_TumorSymlinkToPreProcessingOutputBam:
        input: bam=WGS_preprocessing(expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bam",sample=TUMOR_SAMPLES)),bai=WGS_preprocessing(expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bai",sample=TUMOR_SAMPLES))
        output: bam="out/tumor/variant_calling/{sample}.analysis_ready.bam",bai="out/tumor/variant_calling/{sample}.analysis_ready.bam.bai"
        params: n="1", mem_per_cpu="4", R="'span[hosts=1]'",o="out/logs/variant_calling/symlink.out",eo="out/logs/variant_calling/symlink.err",J="symlink"
        run: 
            in_bam_path = os.path.abspath(str(input.bam))
            in_bai_path = os.path.abspath(str(input.bai))
            command = "ln -s {} {}; ln -s {} {}".format(in_bam_path, output.bam, in_bai_path, output.bai)
            shell(command)
    rule var_00_SymlinkToPreProcessingOutputBam:
        input: bam=WGS_preprocessing(expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bam",sample=NORMAL_SAMPLES)),bai=WGS_preprocessing(expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bai",sample=NORMAL_SAMPLES))
        output: bam=expand("out/{study_group}/variant_calling/{{sample}}.analysis_ready.bam",study_group=['normal'] if config['user_defined_workflow']['genome_personalization_module']['data_is_matched_tumor_normal'] else ['experiment']),bai=expand("out/{study_group}/variant_calling/{{sample}}.analysis_ready.bam.bai",study_group=['normal'] if config['user_defined_workflow']['genome_personalization_module']['data_is_matched_tumor_normal'] else ['experiment'])
        params: n="1", mem_per_cpu="4", R="'span[hosts=1]'",o="out/logs/variant_calling/symlink.out",eo="out/logs/variant_calling/symlink.err",J="symlink"
        run: 
            in_bam_path = os.path.abspath(str(input.bam))
            in_bai_path = os.path.abspath(str(input.bai))
            command = "ln -s {} {}; ln -s {} {}".format(in_bam_path, output.bam, in_bai_path, output.bai)
            shell(command)
elif just_ran_WGS_preprocessing:
    rule var_00_TumorSymlinkToPreProcessingOutputBam:
        input: bam=expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bam",sample=TUMOR_SAMPLES),bai=expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bai",sample=TUMOR_SAMPLES)
        output: bam="out/tumor/variant_calling/{sample}.analysis_ready.bam",bai="out/tumor/variant_calling/{sample}.analysis_ready.bam.bai"
        params: n="1", mem_per_cpu="4", R="'span[hosts=1]'",o="out/logs/variant_calling/symlink.out",eo="out/logs/variant_calling/symlink.err",J="symlink"
        run: 
            in_bam_path = os.path.abspath(str(input.bam))
            in_bai_path = os.path.abspath(str(input.bai))
            command = "ln -s {} {}; ln -s {} {}".format(in_bam_path, output.bam, in_bai_path, output.bai)
            shell(command)
    rule var_00_SymlinkToPreProcessingOutputBam:
        input: bam=expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bam",sample=NORMAL_SAMPLES),bai=expand("out/WGS/{sample}.aligned_sorted_ubam-merged_RG-merged_dedup_fixedtags_BQSR.analysis_ready.bai",sample=NORMAL_SAMPLES)
        output: bam=expand("out/{study_group}/variant_calling/{{sample}}.analysis_ready.bam",study_group=['normal'] if config['user_defined_workflow']['genome_personalization_module']['data_is_matched_tumor_normal'] else ['experiment']),bai=expand("out/{study_group}/variant_calling/{{sample}}.analysis_ready.bam.bai",study_group=['normal'] if config['user_defined_workflow']['genome_personalization_module']['data_is_matched_tumor_normal'] else ['experiment'])
        params: n="1", mem_per_cpu="4", R="'span[hosts=1]'",o="out/logs/variant_calling/symlink.out",eo="out/logs/variant_calling/symlink.err",J="symlink"
        run:
            in_bam_path = os.path.abspath(str(input.bam))
            in_bai_path = os.path.abspath(str(input.bai))
            command = "ln -s {} {}; ln -s {} {}".format(in_bam_path, output.bam, in_bai_path, output.bai)
            shell(command)

#elif 'bam' in input_file_format:
#    print(input_file_format)
else: # input is an analysis-ready BAM
    PROCESSED_BAM_DICT=dict(config['input_files']['genome_personalization_module']['bam_inputs'])
    all_bams_preprocessed = True
    for sample_name in PROCESSED_BAM_DICT.keys():
        if not PROCESSED_BAM_DICT[sample_name]['pre-processing_already_complete']: all_bams_preprocessed = False
    assert('bam' in input_file_format and all_bams_preprocessed), "ERROR: Variant calling is turned on, but WGS/WES pre-processing is turned off. Therefore all input files must be coordinate-sorted, duplicate-marked, analysis-ready BAMs. Please ensure that this is the case, and if so that the corresponding parameters are set (i.e. input_files->genome_personalization_module->input_file_format = 'bam'; input_files->genome_personalization_module->bam_inputs-><sample>->pre-processing_already_complete = true). Please also double check that pre-processing was not disabled erroneously."
    rule var_00_SymlinkToUserPreprocessedBam:
        input: bam=lambda wildcards: os.path.abspath(config['input_files']['genome_personalization_module']['bam_inputs'][wildcards.sample]['bam_file']),bai=lambda wildcards: os.path.abspath(config['input_files']['genome_personalization_module']['bam_inputs'][wildcards.sample]['bai_file'])
        output: bam="out/{study_group}/variant_calling/{sample}.analysis_ready.bam",bai="out/{study_group}/variant_calling/{sample}.analysis_ready.bai"
        params: n="1", mem_per_cpu="4", R="'span[hosts=1]'",o="out/logs/variant_calling/symlink.out",eo="out/logs/variant_calling/symlink.err",J="symlink"
        shell: "ln -s {input.bam} {output.bam}; ln -s {input.bai} {output.bai}"

rule var_germ_01_CallGermlineVariantsPerInterval:
    input: bam="out/{study_group}/variant_calling/{sample}.analysis_ready.bam", interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list"
    output: temp("out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.g.vcf")
    params: n="4", mem_per_cpu="3", R="'span[hosts=1] rusage[mem=3]'", \
        o="out/logs/variant_calling/intervals/vcf_{interval}.out", eo="out/logs/variant_calling/intervals/vcf_{interval}.err", \
        J="generate_vcf_{interval}"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx12g' HaplotypeCaller -R {STOCK_GENOME_FASTA} -I {input.bam} -O {output} -L {input.interval_list} -ERC GVCF"

rule var_germ_02_GenotypeTumorSamplePerInterval:
    input: gvcf="out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.g.vcf", interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list"
    output: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.genotyped.vcf"
    params: n="2", mem_per_cpu="32", R="'span[hosts=1] rusage[mem=32]'", \
        o="out/logs/variant_calling/intervals/genotype_gvcfs_{interval}.out", eo="out/logs/variant_calling/intervals/genotype_gvcfs_{interval}.err", \
        J="genotype_gvcfs"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options '-Xmx64g' GenotypeGVCFs -R {STOCK_GENOME_FASTA} \
              -V {input.gvcf} -L {input.interval_list} -O {output}"

rule var_germ_TMP_MergeHardFilteredIntervalWiseVCFs:
    input: expand("out/{{study_group}}/variant_calling/HTC-scattered/{{sample}}.HTC.{interval}.hardfiltered.vcf.gz",interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])
    output: "out/{study_group}/variant_calling/{sample}.hardfiltered.vcf"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/merge_vcfs.out", eo="out/logs/merge_vcfs.err", \
        J="merge_vcfs"
    conda: "envs/gatk4.yaml"
    shell: "picard MergeVcfs $(echo '{input}' | sed -r 's/[^ ]+/I=&/g') O={output}"

rule var_germ_TMP_merge_genotyped:
    input: vcf=expand("out/{{study_group}}/variant_calling/HTC-scattered/{{sample}}.HTC.{interval}.genotyped.vcf",interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])
    output: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.genotyped.merged.vcf"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/merge_genotyped.out", eo="out/logs/merge_genotyped.err", \
        J="merge_genotyped"
    conda: "envs/gatk4.yaml"
    shell: "picard MergeVcfs $(echo '{input.vcf}' | sed -r 's/[^ ]+/I=&/g') O={output}"

is_whole_genome_or_exome = config['input_files']['genome_personalization_module']['whole_genome_or_exome']
if is_whole_genome_or_exome:
    rule var_germ_TMP_CNN1D_ScoreVariants:
        input: vcf="out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.genotyped.merged.vcf"
        output: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.genotyped.merged.1D_CNN-scored.vcf.gz"
        params: n="16", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/variant_calling/intervals/score_variants.out", eo="out/logs/variant_calling/intervals/score_variants.err", \
            J="score_variants"
        singularity: "docker://broadinstitute/gatk:4.1.4.1"
        shell: "gatk --java-options '-Xmx128g' CNNScoreVariants -R {STOCK_GENOME_FASTA} -V {input.vcf} -O {output}"

    rule var_germ_TMP_CNN2D_ScoreVariants:
        input: bam="out/{study_group}/variant_calling/{sample}.analysis_ready.bam",vcf="out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.genotyped.merged.vcf"
        output: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.genotyped.merged.scored.vcf.gz"
        params: n="16", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/variant_calling/intervals/score_variants.out", eo="out/logs/variant_calling/intervals/score_variants.err", \
            J="score_variants"
        singularity: "docker://broadinstitute/gatk:4.1.4.1"
        shell: "gatk --java-options '-Xmx128g' CNNScoreVariants -R {STOCK_GENOME_FASTA} -I {input.bam} -V {input.vcf} -O {output} --tensor-type read_tensor"
    rule var_germ_03_CNN2D_ScoreVariants:
        input: bam="out/{study_group}/variant_calling/{sample}.analysis_ready.bam",vcf="out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.genotyped.vcf",interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list"
        output: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.CNN-scored.vcf.gz"
        params: n="1", mem_per_cpu="16", R="'span[hosts=1] rusage[mem=16]'", \
            o="out/logs/variant_calling/intervals/score_variants_{interval}.out", eo="out/logs/variant_calling/intervals/score_variants_{interval}.err", \
            J="score_variants_{interval}"
        singularity: "docker://broadinstitute/gatk:4.1.4.1"
        shell: "gatk --java-options '-Xmx8g' CNNScoreVariants -R {STOCK_GENOME_FASTA} -I {input.bam} -V {input.vcf} -O {output} -L {input.interval_list} --tensor-type read_tensor"

    rule var_germ_03_TMP_CNN1D_ScoreVariants:
        input: vcf="out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.genotyped.vcf",interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list"
        output: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.1D_CNN-scored.vcf.gz"
        params: n="1", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/variant_calling/intervals/score_variants_{interval}.out", eo="out/logs/variant_calling/intervals/score_variants_{interval}.err", \
            J="score_variants_{interval}"
        singularity: "docker://broadinstitute/gatk:4.1.4.1"
        shell: "gatk --java-options '-Xmx8g' CNNScoreVariants -R {STOCK_GENOME_FASTA} -V {input.vcf} -O {output} -L {input.interval_list}"

    HAPMAP=config['parameters']['genome_personalization_module']['variant_calling']['resources']['germline']['snps_db']
    MILLS=config['parameters']['genome_personalization_module']['variant_calling']['resources']['germline']['indels_db']
    SNP_TRANCHE=config['parameters']['genome_personalization_module']['variant_calling']['advanced']['snp_tranche']
    INDEL_TRANCHE=config['parameters']['genome_personalization_module']['variant_calling']['advanced']['indel_tranche']
    rule var_germ_04_AssignVariantTranches:
        input: vcf="out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.CNN-scored.vcf.gz",interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list"
        output: temp("out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.CNN-scored.tranched.vcf.gz")
        params: n="1", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", \
            o="out/logs/assign_tranches.out", eo="out/logs/assign_tranches.err", \
            J="assign_tranches"
        singularity: "docker://broadinstitute/gatk:4.1.4.1"
        shell: "gatk --java-options '-Xmx4g' FilterVariantTranches -V {input.vcf} --resource {HAPMAP} --resource {MILLS} --info-key CNN_2D --snp-tranche {SNP_TRANCHE} --indel-tranche {INDEL_TRANCHE} --invalidate-previous-filters -O {output} -L {input.interval_list}"

else:
    
    rule var_germ_0304_HardFilterVariants:
        input: vcf="out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.genotyped.vcf",interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list"
        output: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.hardfiltered.vcf.gz"
        params: n="1", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/variant_calling/intervals/hardfilter_variants_{interval}.out", eo="out/logs/variant_calling/intervals/hardfilter_variants_{interval}.err", \
            J="hardfilter_variants_{interval}"
        conda: "envs/gatk4.yaml"
        shell: "gatk --java-options '-Xmx8g' VariantFiltration -V {input.vcf} -L {input.interval_list} --filter-expression 'QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'HardFiltered' -O {output}" 

"""
rule var_germ_05_FilterNonpassingGermlineVariants:
    input: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.CNN-scored.tranched.vcf.gz"
    output: temp("out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.CNN-scored.tranched.filtered.vcf")
    params: n="1", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", \
        o="out/logs/filter_germline.out", eo="out/logs/filter_germline.err", \
        J="filter_germline"
    conda: "envs/bcftools.yaml"
    shell: "bcftools view -f PASS {input} > {output}"

rule var_germ_06_MergeIntervalWiseVCFs:
    input: expand("out/{{study_group}}/variant_calling/HTC-scattered/{{sample}}.HTC.{interval}.CNN-scored.tranched.filtered.vcf",interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])
    output: "out/{study_group}/variant_calling/{sample}.germline_finished.vcf"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/merge_vcfs.out", eo="out/logs/merge_vcfs.err", \
        J="merge_vcfs"
    conda: "envs/gatk4.yaml"
    shell: "picard MergeVcfs $(echo '{input}' | sed -r 's/[^ ]+/I=&/g') O={output}"
"""

rule var_germ_05_FilterNonpassingGermlineVariants:
    input: "out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.{filter_mode}.vcf.gz"
    output: temp("out/{study_group}/variant_calling/HTC-scattered/{sample}.HTC.{interval}.{filter_mode}.filtered.vcf")
    params: n="1", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", \
        o="out/logs/filter_germline.out", eo="out/logs/filter_germline.err", \
        J="filter_germline"
    conda: "envs/bcftools.yaml"
    shell: "bcftools view -f PASS {input} > {output}"

rule var_germ_06_MergeIntervalWiseVCFs:
    input: expand("out/{{study_group}}/variant_calling/HTC-scattered/{{sample}}.HTC.{interval}.{filter_mode}.filtered.vcf",interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)],filter_mode='CNN-scored.tranched' if is_whole_genome_or_exome else 'hardfiltered')
    output: "out/{study_group}/variant_calling/{sample}.germline_finished.vcf"
    wildcard_constraints: filter_mode="CNN-scored.tranched|hardfiltered"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/merge_vcfs.out", eo="out/logs/merge_vcfs.err", \
        J="merge_vcfs"
    conda: "envs/gatk4.yaml"
    shell: "picard MergeVcfs $(echo '{input}' | sed -r 's/[^ ]+/I=&/g') O={output}"

rule var_germ_06b_ConsolidateSampleNamesForMerge:
    input: "out/{study_group}/variant_calling/{sample}.germline_finished.vcf"
    output: temp("out/{study_group}/variant_calling/{sample}.germline_finished.name_consolidated.vcf.gz"),name_txt=temp("out/{study_group}/variant_calling/{sample}.cohort_name.txt")
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/consolidate_names.out", eo="out/logs/consolidate_names.err", \
        J="consolidate_names",int_vcf="out/{study_group}/variant_calling/{sample}.germline_finished.name_consolidated.vcf"
    conda: "envs/bcftools.yaml"
    shell: "echo '{COHORT}_{wildcards.study_group}' > {output.name_txt}; bcftools reheader -s {output.name_txt} {input} > {params.int_vcf}; bgzip {params.int_vcf}; tabix -p vcf {params.int_vcf}.gz"

rule var_germ_07t_CombineTumorSampleGermlineVCFs:
    input: expand("out/tumor/variant_calling/{sample}.germline_finished.name_consolidated.vcf.gz",sample=TUMOR_SAMPLES)
    output: "out/tumor/variant_calling/{cohort}.tumor.germline_finished.vcf.gz"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/merge_tumor_vcfs.out", eo="out/logs/merge_tumor_vcfs.err", \
        J="merge_tumor_vcfs", int_vcf="out/tumor/variant_calling/{cohort}.tumor.germline_finished.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat -a -D {input} > {params.int_vcf}; bgzip {params.int_vcf}; tabix -p vcf {params.int_vcf}.gz"
    #shell: "picard MergeVcfs $(echo '{input}' | sed -r 's/[^ ]+/I=&/g') O={output}"

rule var_germ_07n_CombineSampleGermlineVCFs:
    input: expand("out/{{study_group}}/variant_calling/{sample}.germline_finished.name_consolidated.vcf.gz",sample=NORMAL_SAMPLES)
    output: "out/{study_group}/variant_calling/{cohort}.{study_group}.germline_finished.vcf.gz","out/{study_group}/variant_calling/{cohort}.{study_group}.germline_finished.vcf.gz.tbi"
    wildcard_constraints: study_group='control|experiment|normal'
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/merge_normal_vcfs.out", eo="out/logs/merge_normal_vcfs.err", \
        J="merge_normal_vcfs", int_vcf="out/{study_group}/variant_calling/{cohort}.{study_group}.germline_finished.vcf"
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat -a -D {input} > {params.int_vcf}; bgzip {params.int_vcf}; tabix -p vcf {params.int_vcf}.gz"


## PINDEL RULES ##

rule CreatePindelConfigFile:
    input: ANALYSIS_READY_BAMFILES
    output: "out/{study_group}/variant_calling/pindel/pindel_config_file.txt"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/create_pindel_configfile.out", eo="out/logs/create_pindel_configfile.err", \
        J="pindel_config"
    run: 
        for sample in ALL_SAMPLES:
            bam_file = "out/{study_group}/variant_calling/{}.analysis_ready.bam".format(sample) if sample in TUMOR_SAMPLES else "out/{study_group}/variant_calling/{}.analysis_ready.bam".format(sample)
            insert_length = config['input_files']['genome_personalization_module'][input_file_format+'_inputs'][sample]['insert_size']
            open(output[0],'a').write("{}\t{}\t{}\n".format(bam_file,insert_length,sample))

try:
    fa_index_file = open('{}.fai'.format(STOCK_GENOME_FASTA))
except FileNotFoundError:
    print("Genome FASTA is not indexed! Please run 'samtools faidx <FASTA>'.")
CHROMOSOMES = [x.split('\t')[0] for x in fa_index_file.readlines()]
rule CallLongIndelsAndSVsWithPindel:
    input: config="out/{study_group}/variant_calling/pindel/pindel_config_file.txt", ref_fasta=STOCK_GENOME_FASTA
    output: "out/{study_group}/variant_calling/pindel/scattered/{cohort}.{chr}.pindel"
    params: n="4", mem_per_cpu="16", R="'span[hosts=1] rusage[mem=16]'", \
        o="out/logs/variant_calling/chr-wise/pindel_{chr}.out", eo="out/logs/variant_calling/chr-wise/pindel_{chr}.err", \
        J="pindel"
    conda: "envs/pindel.yaml"
    shell: "pindel -T {params.n} -f {input.ref_fasta} -i {input.config} -c {wildcards.chr} -o {output}"

rule Pindel2Vcf:
    input: pindel="out/{study_group}/variant_calling/pindel/scattered/{cohort}.{chr}.pindel_TD", ref_fasta=STOCK_GENOME_FASTA
    output: "out/{study_group}/variant_calling/pindel/scattered/{cohort}.{chr}.pindel_TD.vcf"
    params: n="4", mem_per_cpu="16", R="'span[hosts=1] rusage[mem=16]'", \
        o="out/logs/variant_calling/chr-wise/pindel2vcf_{chr}.out", eo="out/logs/variant_calling/chr-wise/pindel2vcf_{chr}.err", \
        J="pindel2vcf"
    conda: "envs/pindel.yaml"
    shell: "pindel2vcf -T {params.n} -p {input.pindel} -r {input.ref_fasta} -R-c {wildcards.chr} -G"

## PINDEL development ongoing ##

## BEGIN Mutect2 RULES ##

SOMATIC_SAMPLE_DICT=dict()
SOMATIC_SAMPLE_DICT['tumor']=TUMOR_SAMPLES
SOMATIC_SAMPLE_DICT['experiment']=NONMATCHED_SAMPLES
SOMATIC_SAMPLE_DICT['normal']=[]
GNOMAD_AF=config['parameters']['genome_personalization_module']['variant_calling']['resources']['somatic']['germline_population_db']
if config['user_defined_workflow']['genome_personalization_module']['data_is_matched_tumor_normal']:
    rule var_som_01_Mutect2_matched_tumor_normal:
        input: tumor=expand("out/tumor/variant_calling/{sample}.analysis_ready.bam",sample=TUMOR_SAMPLES),normal=expand("out/normal/variant_calling/{sample}.analysis_ready.bam",sample=NORMAL_SAMPLES),ref=STOCK_GENOME_FASTA,interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list",gnomad_af=GNOMAD_AF
        output: vcf=temp("out/{study_group}/variant_calling/Mutect2-scattered/{cohort}.mutect2.{interval}.vcf"),stats=temp("out/{study_group}/variant_calling/Mutect2-scattered/{cohort}.mutect2.{interval}.vcf.stats")
        params: n="4", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", \
            o="out/logs/variant_calling/intervals/mutect2_{interval}.out", eo="out/logs/variant_calling/intervals/mutect2_{interval}.err", \
            J="mutect2_matched"
        conda: "envs/gatk4.1.5.0.yaml"
        shell: "gatk --java-options '-Xmx16g' Mutect2 -R {input.ref} -O {output.vcf} \
                  $(echo '{input.tumor}' | sed -r 's/[^ ]+/-I &/g') \
                  $(echo '{input.normal}' | sed -r 's/[^ ]+/-I &/g') \
                  $(echo '{NORMAL_SAMPLES}' | sed -r 's/[^ ]+/-normal &/g') \
                  --germline-resource {input.gnomad_af} \
                  --native-pair-hmm-threads {params.n} \
                  -L {input.interval_list}"
elif 'somatic' in VARIANT_CALLING_MODES:
    rule var_som_01_Mutect2_without_matched_normal:
        input: sample=expand("out/{{study_group}}/variant_calling/{sample}.analysis_ready.bam",sample=NONMATCHED_SAMPLES),ref=STOCK_GENOME_FASTA,interval_list="out/{study_group}/variant_calling/intervals/{interval}-scattered.interval_list",gnomad_af=GNOMAD_AF
        output: vcf="out/{study_group}/variant_calling/Mutect2-scattered/{cohort}.mutect2.{interval}.vcf",stats="out/{study_group}/variant_calling/Mutect2-scattered/{cohort}.mutect2.{interval}.vcf.stats"
        params: n="4", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", \
            o="out/logs/variant_calling/intervals/mutect2_{interval}.out", eo="out/logs/variant_calling/intervals/mutect2_{interval}.err", \
            J="mutect2_unmatched"
        conda: "envs/gatk4.1.5.0.yaml"
        shell: "gatk --java-options '-Xmx16g' Mutect2 -R {input.ref} -O {output.vcf} \
                $(echo '{input.sample}' | sed -r 's/[^ ]+/-I &/g') \
                --germline-resource {input.gnomad_af} -L {input.interval_list}"

rule var_som_02_MergeScatteredMutect2VCFs:
    input: vcf=expand("out/{{study_group}}/variant_calling/Mutect2-scattered/{{cohort}}.mutect2.{interval}.vcf",interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])
    output: "out/{study_group}/variant_calling/{cohort}.mutect2.vcf"
    params: n="1", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/merge_somatic.out", eo="out/logs/merge_somatic.err", \
            J="merge_somatic"
    conda: "envs/gatk4.yaml"
    shell: "picard MergeVcfs -Xmx8g \
              $(echo '{input.vcf}' | sed -r 's/[^ ]+/I=&/g') \
              O={output}"

rule var_som_02a_MergeScatteredMutect2Stats:
    input: stats=expand("out/{{study_group}}/variant_calling/Mutect2-scattered/{{cohort}}.mutect2.{interval}.vcf.stats",interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])
    output: "out/{study_group}/variant_calling/{cohort}.mutect2.vcf.stats"
    params: n="1", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/merge_mutect_stats.out", eo="out/logs/merge_mutect_stats.err", \
            J="merge_mutect_stats"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options -Xmx8g MergeMutectStats -O {output} \
              $(echo '{input.stats}' | sed -r 's/[^ ]+/-stats &/g')"


rule var_som_00a_GetTumorPileupSummaries:
    input: tumor=lambda wildcards: expand("out/{{study_group}}/variant_calling/{sample}.analysis_ready.bam",sample=SOMATIC_SAMPLE_DICT[wildcards.study_group]),gnomad_af=config['parameters']['genome_personalization_module']['variant_calling']['resources']['somatic']['germline_population_db'],intervals=config['parameters']['genome_personalization_module']['variant_calling']['resources']['wgs_calling_regions'],ref_fasta=STOCK_GENOME_FASTA
    output: "out/{study_group}/variant_calling/{cohort}_tumor.pileup_summaries.table"
    params: n="16", mem_per_cpu="2", R="'span[hosts=1] rusage[mem=2]'", \
        o="out/logs/get_pileup_summaries.out", eo="out/logs/get_pileup_summaries.err", \
        J="get_pileup_summaries"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options -Xmx32g GetPileupSummaries -R {input.ref_fasta} \
                $(echo '{input.tumor}' | sed -r 's/[^ ]+/-I &/g') \
                -V {input.gnomad_af} -L {input.intervals} -O {output}"

rule var_som_00b_CalculateContamination:
    input: "out/{study_group}/variant_calling/{cohort}_tumor.pileup_summaries.table"
    output: "out/{study_group}/variant_calling/{cohort}_tumor.contamination.table"
    params: n="4", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", \
        o="out/logs/calculate_contamination.out", eo="out/logs/calculate_contamination.err", \
        J="calculate_contamination"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options -Xmx32g CalculateContamination -I {input} -O {output}"

rule var_som_03_ScoreVariantsWithFilterMutectCalls:
    input: vcf="out/{study_group}/variant_calling/{cohort}.mutect2.vcf",contam_table="out/{study_group}/variant_calling/{cohort}_tumor.contamination.table",stats="out/{study_group}/variant_calling/{cohort}.mutect2.vcf.stats",ref=STOCK_GENOME_FASTA
    output: "out/{study_group}/variant_calling/{cohort}.mutect2.scored.vcf.gz"
    params: n="1", mem_per_cpu="16", R="'span[hosts=1] rusage[mem=16]'", \
        o="out/logs/filter_mutect2.out", eo="out/logs/filter_mutect2.err", \
        J="filter_mutect2"
    conda: "envs/gatk4.yaml"
    shell: "gatk --java-options -Xmx16g FilterMutectCalls -V {input.vcf} -R {input.ref} --contamination-table {input.contam_table} -stats {input.stats} -O {output}" 

rule var_som_04_FilterNonPassingSomaticVariants:
    input: "out/{study_group}/variant_calling/{cohort}.mutect2.scored.vcf.gz"
    output: "out/{study_group}/variant_calling/{cohort}.mutect2.scored.filtered.vcf"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", \
        o="out/logs/filter_somatic.out", eo="out/logs/filter_somatic.err", \
        J="filter_somatic"
    conda: "envs/bcftools.yaml"
    shell: "bcftools view -f PASS {input} > {output}"


rule var_som_05_SeparatePassingSomaticVariantsBySample:
    input: "out/{study_group}/variant_calling/{cohort}.mutect2.scored.filtered.vcf"
    output: "out/{study_group}/variant_calling/{cohort}.{sample}.mutect2.scored.filtered.vcf.gz"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", \
        o="out/logs/variant_calling/somatic_separate_bySample.out", eo="out/logs/variant_calling/somatic_separate_bySample.err", \
        J="somatic_filter+separate", samples=lambda wildcards:SOMATIC_SAMPLE_DICT[wildcards.study_group]
    conda: "envs/gatk4.yaml"
    shell: "for s in {params.samples}; do gatk --java-options -Xmx4g SelectVariants --sample-name $s --variant {input} --output out/{wildcards.study_group}/variant_calling/{COHORT}.$s.mutect2.scored.filtered.vcf.gz; done"

rule var_som_06_ConsolidateSampleNamesForMerge:
    input: "out/{study_group}/variant_calling/{cohort}.{sample}.mutect2.scored.filtered.vcf.gz"
    output: vcf=temp("out/{study_group}/variant_calling/{cohort}.mutect2.scored.filtered.{sample}.name_consolidated.vcf.gz"),name_txt=temp("out/{study_group}/variant_calling/{sample}.cohort_name.{cohort}.txt"),idx=temp("out/{study_group}/variant_calling/{cohort}.mutect2.scored.filtered.{sample}.name_consolidated.vcf.gz.tbi")
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/consolidate_names.out", eo="out/logs/consolidate_names.err", \
        J="consolidate_names",int_vcf="out/{study_group}/variant_calling/{cohort}.mutect2.scored.filtered.{sample}.name_consolidated.vcf"
    conda: "envs/bcftools.yaml"
    shell: "echo '{COHORT}_{wildcards.study_group}' > {output.name_txt}; bcftools reheader -s {output.name_txt} {input} > {output.vcf}; tabix -p vcf {output.vcf}"

rule var_som_07_CombineSomaticVCFs:
    input: vcf=lambda wildcards:expand("out/{{study_group}}/variant_calling/{{cohort}}.mutect2.scored.filtered.{sample}.name_consolidated.vcf.gz",sample=SOMATIC_SAMPLE_DICT[wildcards.study_group]),idx=lambda wildcards: expand("out/{{study_group}}/variant_calling/{{cohort}}.mutect2.scored.filtered.{sample}.name_consolidated.vcf.gz.tbi",sample=SOMATIC_SAMPLE_DICT[wildcards.study_group])
    output: "out/{study_group}/variant_calling/{cohort}.somatic_finished.vcf.gz"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
        o="out/logs/merge_tumor_vcfs.out", eo="out/logs/merge_tumor_vcfs.err", \
        J="merge_tumor_vcfs", int_vcf="out/{study_group}/variant_calling/{cohort}.somatic_finished.vcf"
    wildcard_constraints: study_group='control|experiment|tumor'
    conda: "envs/bcftools.yaml"
    shell: "bcftools concat -a -D {input.vcf} > {params.int_vcf}; bgzip {params.int_vcf}; tabix -p vcf {params.int_vcf}.gz"

rule var_z_MergeFinishedTumorVCFs:
    input: ["out/tumor/variant_calling/{cohort}.somatic_finished.vcf.gz","out/tumor/variant_calling/{cohort}.tumor.germline_finished.vcf.gz"]
    output: "out/tumor/variant_calling/{cohort}.tumor.variant_calling_finished.vcf.gz"
    params: n="1", mem_per_cpu="8", R="'span[hosts=1] rusage[mem=8]'", \
            o="out/logs/merge_finished_vcfs.out", eo="out/logs/merge_finished_vcfs.err", \
            J="merge_finished_vcfs"
    conda: "envs/gatk4.yaml"
    shell: "picard MergeVcfs $(echo '{input}' | sed -r 's/[^ ]+/I=&/g') O={output}"

rule var_z_FinishNonTumorVCF:
    input: ["out/{study_group}/variant_calling/{cohort}.{study_group}.germline_finished.vcf.gz","out/{study_group}/variant_calling/{cohort}.somatic_finished.vcf.gz"] if 'somatic' in VARIANT_CALLING_MODES and TUMOR_SAMPLES==[] else ["out/{study_group}/variant_calling/{cohort}.{study_group}.germline_finished.vcf.gz"]
    #input: vcf="out/{study_group}/variant_calling/{cohort}.{study_group}.germline_finished.vcf.gz",tbi="out/{study_group}/variant_calling/{cohort}.{study_group}.germline_finished.vcf.gz.tbi"
    output: vcf="out/{study_group}/variant_calling/{cohort}.{study_group}.variant_calling_finished.vcf.gz"#,tbi="out/{study_group}/variant_calling/{cohort}.{study_group}.variant_calling_finished.vcf.gz.tbi"
    params: n="1", mem_per_cpu="4", R="'span[hosts=1]'", \
            o="out/logs/merge_finished_vcfs.out", eo="out/logs/merge_finished_vcfs.err", \
            J="merge_finished_vcfs"
    wildcard_constraints: study_group='control|experiment|normal'
    conda: "envs/gatk4.yaml"
    shell: "picard MergeVcfs $(echo '{input}' | sed -r 's/[^ ]+/I=&/g') O={output}"
