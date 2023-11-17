## Setup

- Run on **terra**.

    ```plain
    export PATH=/opt/common/bic/miniconda3/bin:$PATH
    source activate /home/wilson/.conda/envs/snakemakev5 ## Manda already set up snakemake env successfully
    module load singularity/3.7.1

    #conda deactivate
    ```


- Nick copied full paths to example input & output from lilac to juno in:

    ```
    /juno/res/mix/Cache/2023-09-14/2023_10_test_data
    ```
- Stock references and resources are copied (full paths) from lilac to juno in:

    ```
    /juno/depot/custom/PG2
    ```


- Copy original sample config file to project directory. Keep original for reference.
    ```plain
    cp /juno/res/mix/Cache/2023-09-14/2023_10_test_data/data/kentsis/ProteomeGenerator2/configfiles/2023-Oct-12-KasuminonDAC.yaml \
        /juno/work/bic/byrne/pg2/20231023_new_cmd/2023-Oct-12-KasuminonDAC.yaml.ORIGINAL 
    ```

- Make another copy and modify the new one to point to the correct input files.
    ```plain
    cp /juno/work/bic/byrne/pg2/20231023_new_cmd/2023-Oct-12-KasuminonDAC.yaml.ORIGINAL \
        /juno/work/bic/byrne/pg2/20231023_new_cmd/2023-Oct-12-KasuminonDAC.yaml

    # Modify /juno/work/bic/byrne/pg2/20231023_new_cmd/2023-Oct-12-KasuminonDAC.yaml
    # NOTE: Since Nick copied full paths to data, in our config some original paths
    # to stock references are prepended with /juno/depot/custom/PG2 and original paths
    # to input (fastqs) are prepended with /juno/res/mix/Cache/2023-09-14/2023_10_test_data
    ```
    **NOTE:** both original and modified configs are also in [ProteomeGenerator2/BIC_config](../BIC_config)


- Make a project-specific copy of `run.sh` and modify it to point to your project directory, config file, cluster job submission command, target file and log file. 
    ```plain
    cd /juno/work/bic/byrne/ProteomeGenerator2
    cp run.sh run_KasuminonDAC_test.sh

    #open and modify variables in run_KasuminonDAC_test.sh

    WD=/juno/work/bic/byrne/pg2/20231023_new_cmd                                    
    CONFIG=${WD}/2023-Oct-12-KasuminonDAC.yaml                                      
    TARGET=out/experiment/novel_analysis/proteome_blast.outfmt6                     
    LOG=${WD}/2023-Oct-12-KasuminonDAC.snakemake.out                                
    CLUSTER="bsub -J {params.J} -n {params.n} -R 'span[hosts=1] rusage[mem={params.mem_per_cpu}]' -W 144:00 -o {params.o} -eo {params.eo}"
    
    ...
    ```

- For our test run, the final command was:
    ```plain
    ProteomeGenerator2> snakemake --snakefile ProteomeGenerator2.py \
       --cluster "bsub -J {params.J} -n {params.n} -R 'span[hosts=1] rusage[mem={params.mem_per_cpu}]' -W 144:00 -o {params.o} -eo {params.eo}" \
       -j 100 \
       -k \
       --ri \
       --latency-wait 30 \
       --configfile /juno/work/bic/byrne/pg2/20231023_new_cmd/2023-Oct-12-KasuminonDAC.yaml \
       --use-conda \
       --use-singularity \
       --singularity-args "--bind /juno:/juno,/work:/work,/home:/home,/scratch:/scratch" \
       out/experiment/novel_analysis/proteome_blast.outfmt6 \
      > /juno/work/bic/byrne/pg2/20231023_new_cmd/2023-Oct-12-KasuminonDAC.snakemake.out 2>&1
    ```
---

## Run 
```plain
# move to repo directory if not already there
cd /juno/work/bic/byrne/ProteomeGenerator2

# dry-run
./run_KasuminonDAC_test.sh -n --quiet

# run in background as entire pipeline will take a few days
nohup ./run_KasuminonDAC_test.sh &
```
---

## NOTES -- Changes made from original instructions
- reordered channels in envs/crossmap.yaml
- replaced channel 'bioconda' with 'bioconda/label/cf201901' in envs/transdecoder.yaml
- ~~load perl prior to running (`module load perl/perl-5.22.0`)~~ added channel 'conda-forge' and dependency perl-5.22.0 in envs/transdecoder.yaml
- added /work to singularity bind paths
- put snakemake command in shell script: run.sh; variables will need to be configured in this script but then pipeline can be started and restarted as necessary with simple command `./run.sh`
- increased memory request for job var_germ_03_CNN2D_ScoreVariants 
