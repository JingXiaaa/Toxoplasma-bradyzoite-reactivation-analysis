#!/usr/bin/env bash

BASE="./Cellranger"
GTF="./reference/Tgondii_ME49.gtf"

for sample in EW_0h EW_4h EW_non_invasive EW_8h EW_24h; do
    velocyto run10x \
        "$BASE/count_${sample}/" \
        "$GTF" \
        > "velocyto_${sample}.log" 2>&1 &
done

wait
