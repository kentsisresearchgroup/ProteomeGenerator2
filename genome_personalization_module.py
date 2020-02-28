WD = config['directories']['working_and_output_dir']
workdir: WD
PG2_HOME=config['directories']['PG2_installation_dir']
TMP_DIR=config['directories']['optional']['temp_dir']


STOCK_GENOME_FASTA=config['stock_references']['genome']['fasta']
STOCK_GENOME_GTF=config['stock_references']['genome']['gtf']


COHORT=config['input_files']['genome_personalization_module']['cohort_or_organism_name']

input_file_format = config['input_files']['genome_personalization_module']['input_file_format']
if input_file_format == 'vcf':
    INPUT_VCFS = config['input_files']['genome_personalization_module']['vcf_inputs'].keys()
    INPUT_VCF_GZ_FILES = []
    ALL_INPUT_VCF_SAMPLES = []
    input_vcf_aligned_fastas = []
    INPUT_VCF_GZ_FILES = [config['input_files']['genome_personalization_module']['vcf_inputs'][x]['vcf.gz_file'] for x in INPUT_VCFS]
    for vcf in INPUT_VCFS:
        for name in config['input_files']['genome_personalization_module']['vcf_inputs'][vcf]['sample_names']:
            ALL_INPUT_VCF_SAMPLES.append(name)
    input_vcf_aligned_fastas = [config['input_files']['genome_personalization_module']['vcf_inputs'][x]['aligned_genome']['fasta'] for x in INPUT_VCFS]
    
    for fa in input_vcf_aligned_fastas:
        assert(fa == STOCK_GENOME_FASTA), "Error: At least one of your input vcfs is aligned to a different reference genome than the provided stock ref genome.\n\
            Please use UCSC liftOver or CrossMap to realign each vcf's genome coordinates with the provided stock ref genome.\n\
            In the future, this will be performed automatically."


calling_variants = config['user_defined_workflow']['genome_personalization_module']['variant_calling_submodule']['call_germline_variants_with_GATK4_HaplotypeCaller'] or config['user_defined_workflow']['genome_personalization_module']['variant_calling_submodule']['if_inputs_are_matched_tumor-normal_samples']['call_somatic_variants_with_GATK4_Mutect2']
#calling_variants = False if input_file_format == 'vcf' else True # True if input is 'fastq' or 'bam'

if not calling_variants: assert input_file_format == 'vcf'
if calling_variants:
    TUMOR_SAMPLES=[]
    #NORMAL_SAMPLES=[]
    sample_dict = dict(config['input_files']['genome_personalization_module'][input_file_format+'_inputs'])
    for sample_name in sample_dict.keys():
        if sample_dict[sample_name]['matched_sample_params']['is_matched_sample'] == False or sample_dict[sample_name]['matched_sample_params']['tumor_or_normal']=='tumor': TUMOR_SAMPLES.append(sample_name)
        #elif sample_dict[sample_name]['matched_sample_params']['tumor_or_normal']=='normal':
        #    NORMAL_SAMPLES.append(sample_name)
        #else: assert False, "This should never happen!"

try:
    fa_index_file = open('{}.fai'.format(STOCK_GENOME_FASTA))
except FileNotFoundError:
    print("Genome FASTA is not indexed! Please run 'samtools faidx <FASTA>'.")
CHROMOSOMES = [x.split('\t')[0] for x in fa_index_file.readlines()]
HAPLOTYPES = [1,2] if config['parameters']['genome_personalization_module']['variant_calling']['make_customRef_diploid'] else [1]

rule all:
    input: expand("out/custom_ref/{cohort}_H{htype}.fa",cohort=COHORT,htype=HAPLOTYPES), expand("out/custom_ref/{cohort}_H{htype}.chain",cohort=COHORT,htype=HAPLOTYPES), expand("out/custom_ref/{cohort}_H{htype}.gtf",cohort=COHORT,htype=HAPLOTYPES)


snakemake.utils.makedirs('out/logs/chr-wise')
snakemake.utils.makedirs('out/logs/custom_ref')


subworkflow WGS_variant_calling:
    snakefile: "WGS_variant_calling.py"
    configfile: workflow.overwrite_configfile
    workdir: WD

# For each chromosome, patch the stock reference with variants contained in the (phased) VCF.
# *we use the included dev branch of bcftools because the latest release (1.9) has a bug in bcftools consensus*
BCFTOOLS_DEV = os.path.join(PG2_HOME, config['non-conda_packages']['bcftools-dev'])
if calling_variants:
    rule genome_01_CreateChrWiseCustomRef:
        input: vcf=WGS_variant_calling("out/WGS/variant_calling/cohort/{cohort}.variant_calling_finished.vcf.gz")
        output: fasta=temp("out/custom_ref/chr_split/{cohort}.h-{htype}^{chr}.fa"),chain=temp("out/custom_ref/chr_split/{cohort}.h-{htype}^{chr}.chain")
        params: n="1", R="'rusage[mem=4]'", J="chr-wise_customRef", o="out/logs/chr-wise/{htype}^{chr}.out", eo="out/logs/chr-wise/{htype}^{chr}.err", \
                samples=TUMOR_SAMPLES[0]
        conda: "{PG2_HOME}/envs/bcftools.yaml"
        shell: "samtools faidx {STOCK_GENOME_FASTA} {wildcards.chr} | bcftools consensus -s {params.samples} -H {wildcards.htype} -p {wildcards.htype}_ -c {output.chain} {input.vcf} > {output.fasta}"
else:
    rule genome_00_MergeAllInputVCFs:
        input: vcfs=INPUT_VCF_GZ_FILES
        output: merged_vcf="out/WGS/{cohort}.input_vcfs_merged.vcf.gz"
        conda: "{PG2_HOME}/envs/myenv.yaml"
        params: n="1", R="'rusage[mem=6]'", J="merge_input_vcfs", o="out/logs/merge_input_vcfs.out", eo="out/logs/merge_input_vcfs.err"
        conda: "{PG2_HOME}/envs/bcftools.yaml"
        shell: "picard MergeVcfs -Xmx6g \
            $(echo {input.vcfs} | sed -r 's/[^ ]+/INPUT=&/g') \
            OUTPUT={output.merged_vcf}"
    rule genome_01_CreateChrWiseCustomRef:
        input: merged_vcf="out/WGS/{cohort}.input_vcfs_merged.vcf.gz"
        output: fasta=temp("out/custom_ref/chr_split/{cohort}.h-{htype}^{chr}.fa"),chain=temp("out/custom_ref/chr_split/{cohort}.h-{htype}^{chr}.chain")
        conda: "{PG2_HOME}/envs/bcftools.yaml"
        params: n="1", R="'rusage[mem=4]'", J="chr-wise_customRef", o="out/logs/chr-wise/{htype}^{chr}.out", eo="out/logs/chr-wise/{htype}^{chr}.err",samples=ALL_INPUT_VCF_SAMPLES[1]
        shell: "samtools faidx {STOCK_GENOME_FASTA} {wildcards.chr} | bcftools consensus -s {params.samples} -H {wildcards.htype} -p {wildcards.htype}_ -c {output.chain} {input.merged_vcf} > {output.fasta}"

rule genome_02a_MergeChrWiseCustomRef:
    input: expand("out/custom_ref/chr_split/{{cohort}}.h-{{htype}}^{chr}.fa",chr=CHROMOSOMES),expand("out/custom_ref/chr_split/{{cohort}}.h-{{htype}}^{chr}.chain",chr=CHROMOSOMES)
    output: fasta="out/custom_ref/{cohort}_H{htype}.fa",chain="out/custom_ref/{cohort}_H{htype}.chain"
    params: n="1", R="'rusage[mem=4]'", J="merge_customRef", o="out/logs/merge_customRef.out", eo="out/logs/merge_customRef.err"
    shell: "awk 1 out/custom_ref/chr_split/{wildcards.cohort}.h-{wildcards.htype}^*.fa > {output.fasta}; awk 1 out/custom_ref/chr_split/{wildcards.cohort}.h-{wildcards.htype}^*.chain > {output.chain}"

# Adjust ref annotation coords to account for indels
rule genome_02b_LiftoverAnnotationGTF:
    input: chain="out/custom_ref/{cohort}_H{htype}.chain", vcf_refGtf=STOCK_GENOME_GTF
    output: "out/custom_ref/{cohort}_H{htype}_temp.gtf"
    params: n="1", R="'rusage[mem=4]'", J="LiftoverGTF", o="out/logs/liftover.out", eo="out/logs/liftover.err"
    conda: "{PG2_HOME}/envs/crossmap.yaml"
    shell: "CrossMap.py gff {input.chain} {input.vcf_refGtf} | awk '\{print \"{wildcards.htype}\" $0\}' > {output}"

rule CreateRefSequenceIndex:
    input: STOCK_GENOME_FASTA
    output: STOCK_GENOME_FASTA+".fai"
    params: n="1", R="'span[hosts=1] rusage[mem=18]'", o="out/logs/create_refIdx.out", eo="out/logs/create_refDict.err", J="create_refIdx"
    conda: "{PG2_HOME}/envs/bcftools.yaml"
    shell: "samtools faidx {input}"

rule CreateRefSequenceDict:
    input: STOCK_GENOME_FASTA
    output: STOCK_GENOME_FASTA.strip('fa')+'dict'
    params: n="1", R="'span[hosts=1] rusage[mem=18]'", o="out/logs/create_refDict.out", eo="out/logs/create_refDict.err", J="create_refDict"
    conda: "{PG2_HOME}/envs/bcftools.yaml"
    shell: "picard -Xmx16g CreateSequenceDictionary R={input} O={output}"
