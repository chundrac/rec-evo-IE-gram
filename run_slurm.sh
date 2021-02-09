#!/bin/sh

module use /sapps/etc/modules/start/
module load generic

for i in {1..65}
do
    sbatch run.sh Rscript rate_inference.R $i
done
