#!/usr/bin/sh

#
# Modify and uncomment ALL variables below. Script will not run without all values set!
#

#WD=/juno/work/bic/byrne/pg2/fork_test
#CONFIG=${WD}/2023-Oct-12-KasuminonDAC.yaml
#LOG=${WD}/2023-Oct-12-KasuminonDAC.snakemake.out                                

# vars below will not generally need to be changed after initial configuration
#TARGET=out/experiment/novel_analysis/proteome_blast.outfmt6
#CLUSTER="bsub -J {params.J} -n {params.n} -R 'span[hosts=1] rusage[mem={params.mem_per_cpu}]' -W 144:00 -o {params.o} -eo {params.eo}"
#SINGULARITY="--bind /juno:/juno,/work:/work,/home:/home,/scratch:/scratch"

snakemake --snakefile ProteomeGenerator2.py \
  --configfile "${CONFIG}" \
  --cluster "${CLUSTER}" \
  -j 100 \
  -k \
  --ri \
  --latency-wait 30 \
  --use-conda \
  --use-singularity \
  --singularity-args "${SINGULARITY}" \
  "${TARGET}" \
  "$@" \
  > ${LOG} 2>&1

