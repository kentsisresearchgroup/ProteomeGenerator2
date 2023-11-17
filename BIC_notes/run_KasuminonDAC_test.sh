#!/usr/bin/sh

WD=/juno/work/bic/byrne/pg2/20231023_new_cmd
CONFIG=${WD}/2023-Oct-12-KasuminonDAC.yaml
TARGET=out/experiment/novel_analysis/proteome_blast.outfmt6
LOG=${WD}/2023-Oct-12-KasuminonDAC.snakemake.out
CLUSTER="bsub -J {params.J} -n {params.n} -R 'span[hosts=1] rusage[mem={params.mem_per_cpu}]' -W 144:00 -o {params.o} -eo {params.eo}"

snakemake --snakefile ProteomeGenerator2.py \
  --configfile ${CONFIG} \
  --cluster ${CLUSTER} \
  -j 100 \
  -k \
  --ri \
  --latency-wait 30 \
  --use-conda \
  --use-singularity \
  --singularity-args "--bind /juno:/juno,/work:/work,/home:/home,/scratch:/scratch" \
  "$@" \
  ${TARGET} \
  > ${LOG} 2>&1

