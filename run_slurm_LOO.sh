#!/bin/sh

module use /sapps/etc/modules/start/
module load generic

for i in {1..65}
  for j in {1..148}
    do
      sbatch run.sh Rscript LOO.R $i $j
    done
