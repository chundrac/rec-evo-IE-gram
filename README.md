# rec-evo-IE-gram

Code for paper "Reconstructing the evolution of Indo-European Grammar"

For main analyses, the following script can be run on a SLURM workload manager:

./run_slurm.sh

To generate tree graphics on a SLURM system:

for i in {1..65}
do
    sbatch run.sh Rscript anc_rec.R $i
done

To analyze results:

Rscript summarize_results.R
Rscript summarize_reconstructions_for_comparison.R

All remaining scripts can be run after these steps.
