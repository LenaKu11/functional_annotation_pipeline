# Functional annotation pipeline for genomes on the Life Science Compute Cluster.

Annotation directory, path to genomic sequence and name has to be altered in the file 'config.yaml'
The file 'snakefile' contains the snake-rules, which can be adapted, especially for interproscan (monthly mirror may need to be adapted to the most recent version on the cluster), flags, or versions of programs.

This code has to be run in a tmux session, see https://tmuxcheatsheet.com/ for further information.

To start a session, type:

`tmux new`

And then

`module load miniconda`
`module load snakemake`


Command for a dry mode:

snakemake -npr -j 100 --use-conda --cluster-config cluster.yaml --cluster "sbatch -t {cluster.time} --cpus-per-task {cluster.cpus-per-task} --mem {cluster.mem} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}"  &> dry_run_Cmulti_DATE.out

Main run:

snakemake -j 100 --rerun-incomplete --latency-wait 30 --use-conda --cluster-config cluster.yaml --cluster "sbatch -t {cluster.time} --cpus-per-task {cluster.cpus-per-task} --mem {cluster.mem} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}"  &> main_run_Cmulti_DATE.out


Programs and versions used and tested:



