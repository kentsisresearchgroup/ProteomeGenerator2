CHROMOSOMES=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']

PROTEIN_FA="/juno/depot/custom/PG2/data/kentsis/indexes/GRCh38/gencode.v31.pc_translations.fa"
MUTATION_TYPES=['missense','insertions','deletions','frameshifts']

continuing_after_variant_calling=config['user_defined_workflow']['genome_personalization_module']['variant_calling_submodule']['continuation']['just_finished_variant_calling']

subworkflow WGS_variant_calling:
    snakefile: "WGS_variant_calling.py"
    configfile: workflow.overwrite_configfile
    workdir: WD

snakemake.utils.makedirs('out/logs/novel_analysis')
[snakemake.utils.makedirs('out/{}/novel_analysis/{}'.format(study_group,mutation_type)) for study_group in STUDY_GROUPS for mutation_type in MUTATION_TYPES]

rule IdentifyNovelPeptides:
    input: mq_peps="out/{study_group}/MaxQuant/combined/txt/peptides.txt", ref_db=STOCK_PROTEOME_FASTA
    output: novel_peps="out/{study_group}/novel_analysis/novel_peptides.txt",novelpep_transcript_map="out/{study_group}/novel_analysis/novelPeptide_transcript_map.txt"
    params: n="1", mem_per_cpu="8", R="'rusage[mem=8]'", J="filter_novel_peps", o="out/logs/novel_analysis/{study_group}.filter_novel_peps.out", eo="out/logs/novel_analysis/{study_group}.filter_novel_peps.err"
    shell: "python3 {PG2_HOME}/scripts/generate_novels.py {input.mq_peps} {input.ref_db} {output.novel_peps} {output.novelpep_transcript_map}"

rule AdjustProteomeFastaHeaders:
    input: "out/{study_group}/combined.proteome.unique.fasta"
    output: "out/{study_group}/combined.proteome.unique.headers_adjusted.fasta"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="adjust_proteome_headers", o="out/logs/novel_analysis/{study_group}.adjust_proteome_headers.out", eo="out/logs/novel_analysis/{study_group}.adjust_proteome_headers.err"
    shell: "python3 {PG2_HOME}/scripts/adjust_proteome_fasta_headers.py {input} > {output}"

rule BlastProteome:
    input: proteome="out/{study_group}/combined.proteome.unique.headers_adjusted.fasta", ref_db="/juno/depot/custom/PG2/data/kentsis/indexes/GRCh38/gencode.v31.pc_translations.fa"
    output: "out/{study_group}/novel_analysis/proteome_blast.outfmt6"
    params: n="24", mem_per_cpu="3", R="'span[hosts=1] rusage[mem=3]'", J="blast_proteome", o="out/logs/blast_proteome.out", eo="out/logs/blast_proteome.err"
    conda: "envs/blast.yaml"
    shell: "blastp -num_threads {params.n} -query {input.proteome} -db {input.ref_db} -outfmt 6 -max_target_seqs 5 -evalue 1e-40 > {output}"

#SNPEFF_DB=os.path.splitext(os.path.basename(STOCK_GENOME_GTF))[0]
SNPEFF_DB=os.path.splitext(os.path.basename(PROTEIN_FA))[0]
rule CreateSnpEffDatabase:
    input: ref_gtf=STOCK_GENOME_GTF, genome_fa=STOCK_GENOME_FASTA, protein_fa="/juno/depot/custom/PG2/data/kentsis/indexes/GRCh38/gencode.v31.pc_translations.fa", snpEff_config=os.path.join(PG2_HOME,'utils/snpEff.config')
    output: genome_fa=expand("{pg2_home}/utils/snpEff_dbs/{db_name}/sequences.fa",pg2_home=PG2_HOME,db_name=SNPEFF_DB),gtf=expand("{pg2_home}/utils/snpEff_dbs/{db_name}/genes.gtf",pg2_home=PG2_HOME,db_name=SNPEFF_DB),protein_fa=expand("{pg2_home}/utils/snpEff_dbs/{db_name}/protein.fa",pg2_home=PG2_HOME,db_name=SNPEFF_DB)
    params: n="8", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", J="snpeff_db", o="out/logs/snpeff_db.out", eo="out/logs/snpeff_db.err", db_name=SNPEFF_DB
    conda: "envs/biopython.yaml"
    shell: "ln -s /juno/depot/custom/PG2/data/kentsis/indexes/GRCh38/gencode.v31.pc_translations.fa {PG2_HOME}/utils/snpEff_dbs/{params.db_name}/protein.fa; ln -s {input.ref_gtf} {output.gtf}; ln -s {input.genome_fa} {output.genome_fa}; cat <( echo '{params.db_name}.genome : Human') {input.snpEff_config} > {input.snpEff_config}.new; mv {input.snpEff_config}.new {input.snpEff_config}; snpEff -Xmx32g build -c {input.snpEff_config} -gtf22 -v {params.db_name}"

rule AnnotatePredictedVariantEffects:
    input: vcf="out/{study_group}/variant_calling/{cohort}.{study_group}.variant_calling_finished.vcf.gz" if continuing_after_variant_calling else WGS_variant_calling("out/{study_group}/variant_calling/{cohort}.{study_group}.variant_calling_finished.vcf.gz"),snpEff_config=os.path.join(PG2_HOME,'utils/snpEff.config'), fa=expand("{pg2_home}/utils/snpEff_dbs/{db_name}/sequences.fa",pg2_home=PG2_HOME,db_name=SNPEFF_DB), gtf=expand("{pg2_home}/utils/snpEff_dbs/{db_name}/genes.gtf",pg2_home=PG2_HOME,db_name=SNPEFF_DB), protein_fa=expand("{pg2_home}/utils/snpEff_dbs/{db_name}/protein.fa",pg2_home=PG2_HOME,db_name=SNPEFF_DB)
    output: vcf="out/{study_group}/variant_calling/{cohort}.{study_group}.variant_calling_finished.snpEff.vcf.gz",tbi="out/{study_group}/variant_calling/{cohort}.{study_group}.variant_calling_finished.snpEff.vcf.gz.tbi"
    params: n="16", mem_per_cpu="4", R="'span[hosts=1] rusage[mem=4]'", J="snpEff", o="out/logs/annotate_variants.out", eo="out/logs/annotate_variants.err", \
            db_name=SNPEFF_DB, int_vcf="out/{study_group}/variant_calling/{cohort}.{study_group}.variant_calling_finished.snpEff.vcf"
    conda: "envs/biopython.yaml"
    shell: "snpEff -Xmx64g -c {input.snpEff_config} {params.db_name} {input.vcf} > {params.int_vcf}; bgzip {params.int_vcf}; tabix -p vcf {output.vcf}"

rule SubsetAnnotatedVariantsByChr:
    input: snpEff=expand("out/{{study_group}}/variant_calling/{cohort}.{{study_group}}.variant_calling_finished.snpEff.vcf.gz",cohort=COHORT)
    output: 'out/{study_group}/novel_analysis/snpEff/{chr}.snpEff.vcf'
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="subset_annotatedVCF", o="out/logs/novel_analysis/{study_group}.subset_snpeff.out", eo="out/logs/novel_analysis/{study_group}.subset_snpeff.err"
    conda: "envs/bcftools.yaml"
    shell: "bcftools view {input.snpEff} {wildcards.chr} > {output}"

"""
rule MapMutations_noMS:
    input: annotated_vcf=lambda wildcards:expand("out/{study_group}/novel_analysis/snpEff/{{chr}}.snpEff.vcf",study_group=wildcards.study_group if matched_tumor_normal_data else 'experiment'),proteome="out/{study_group}/combined.proteome.unique.headers_adjusted.fasta",proteome_blast="out/{study_group}/novel_analysis/proteome_blast.outfmt6",ref_db="/juno/depot/custom/PG2/data/kentsis/indexes/GRCh38/gencode.v31.pc_translations.fa"
    output: analysis='out/{study_group}/novel_analysis/{mutation_type}/{chr}.{mutation_type}.analysis'
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="map_{mutation_type}", o="out/logs/novel_analysis/{study_group}.{mutation_type}.out", eo="out/logs/novel_analysis/{study_group}.{mutation_type}.err"
    conda: "envs/biopython.yaml"
    shell: "python {PG2_HOME}/scripts/{wildcards.mutation_type}_noMS_5-10.py {input.proteome} {input.annotated_vcf} {input.proteome_blast} {input.ref_db} > {output.analysis}"
"""


rule MapMutations:
    input: novel_peps="out/{study_group}/novel_analysis/novel_peptides.txt",annotated_vcf=lambda wildcards:expand("out/{study_group}/novel_analysis/snpEff/{{chr}}.snpEff.vcf",study_group=wildcards.study_group if matched_tumor_normal_data else 'experiment'),proteome="out/{study_group}/combined.proteome.unique.headers_adjusted.fasta",proteome_blast="out/{study_group}/novel_analysis/proteome_blast.outfmt6",ref_db="/juno/depot/custom/PG2/data/kentsis/indexes/GRCh38/gencode.v31.pc_translations.fa",novelpep_transcript_map="out/{study_group}/novel_analysis/novelPeptide_transcript_map.txt"
    output: analysis='out/{study_group}/novel_analysis/{mutation_type}/{chr}.{mutation_type}.analysis'
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="map_{mutation_type}", o="out/logs/novel_analysis/{study_group}.{mutation_type}.out", eo="out/logs/novel_analysis/{study_group}.{mutation_type}.err", \
            mutation_mstrg_map='out/{study_group}/novel_analysis/{mutation_type}/{chr}.{mutation_type}.map',mutation_MQevidence_map="out/{study_group}/novel_analysis/{mutation_type}/{chr}.{mutation_type}_MQevidence.map",novelpep_mutation_map="out/{study_group}/novel_analysis/{mutation_type}/{chr}.novelPep_{mutation_type}.map"
    conda: "envs/biopython.yaml"
    shell: "python {PG2_HOME}/scripts/{wildcards.mutation_type}_snpEff.py {input.proteome} {input.annotated_vcf} {input.proteome_blast} {input.ref_db} {input.novelpep_transcript_map} {params.mutation_mstrg_map} {params.mutation_MQevidence_map} {params.novelpep_mutation_map} > {output.analysis}"
    #shell: "python3 {PG2_HOME}/scripts/{wildcards.mutation_type}_snpEff.py {input.proteome} {input.annotated_vcf} {input.proteome_blast} {input.ref_db} {input.novelpep_transcript_map} {params.mutation_mstrg_map} {params.mutation_MQevidence_map} {params.novelpep_mutation_map} > {output.analysis}"


rule AggregateMutations:
    input: expand("out/{{study_group}}/novel_analysis/{{mutation_type}}/{chr}.{{mutation_type}}.analysis", chr=CHROMOSOMES)
    output: mutation_mstrg="out/{study_group}/novel_analysis/{mutation_type}/combined.{mutation_type}.map",mutation_MQevidence="out/{study_group}/novel_analysis/{mutation_type}/combined.{mutation_type}_MQevidence.map",novelPep_mutation="out/{study_group}/novel_analysis/{mutation_type}/combined.novelPep_{mutation_type}.map"
    params: n="1", mem_per_cpu="4", R="'rusage[mem=4]'", J="aggregate_{mutation_type}", o="out/logs/novel_analysis/{study_group}.aggregate_{mutation_type}.out", eo="out/logs/novel_analysis/{study_group}.aggregate_{mutation_type}.err", \
            mutation_mstrg=lambda wildcards:expand('out/{study_group}/novel_analysis/{mutation_type}/{chr}.{mutation_type}.map',chr=[re.search(r'chr[0-9X]+',x).group() for x in [f for f in os.listdir("out/{}/novel_analysis/{}".format(wildcards.study_group,wildcards.mutation_type)) if '.{}.map'.format(wildcards.mutation_type) in f and 'combined' not in f]],study_group=wildcards.study_group,mutation_type=wildcards.mutation_type), \
            mutation_MQevidence=lambda wildcards: expand("out/{study_group}/novel_analysis/{mutation_type}/{chr}.{mutation_type}_MQevidence.map",chr=[re.search(r'chr[0-9X]+',x).group() for x in [f for f in os.listdir("out/{}/novel_analysis/{}".format(wildcards.study_group,wildcards.mutation_type)) if 'MQevidence' in f and 'combined' not in f]],study_group=wildcards.study_group,mutation_type=wildcards.mutation_type), \
            novelpep_mutation=lambda wildcards: expand("out/{study_group}/novel_analysis/{mutation_type}/{chr}.novelPep_{mutation_type}.map",chr=[re.search(r'chr[0-9X]+',x).group() for x in [f for f in os.listdir("out/{}/novel_analysis/{}".format(wildcards.study_group,wildcards.mutation_type)) if 'novelPep' in f and 'combined' not in f]],study_group=wildcards.study_group,mutation_type=wildcards.mutation_type)
    shell: "python3 {PG2_HOME}/scripts/aggregate_mutations.py {params.mutation_mstrg} {params.mutation_MQevidence} {params.novelpep_mutation}"

