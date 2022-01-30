# Functional annotation pipeline on the Life Science Compute Cluster (Vienna).

Annotation directory, path to genomic sequence and name has to be altered in the file 'config.yaml'
The file 'snakefile' contains the snake-rules, which can be adapted, especially for interproscan (monthly mirror may need to be adapted to the most recent version on the cluster), flags, or versions of programs.

This code has to be run in a tmux session, see https://tmuxcheatsheet.com/ for further information.

To start a session, type:

`tmux new`

And then

`module load miniconda`

`module load snakemake`


** Command for a dry run: **

`snakemake -npr -j 100 --use-conda --cluster-config cluster.yaml --cluster "sbatch -t {cluster.time} --cpus-per-task {cluster.cpus-per-task} --mem {cluster.mem} --mail-type {cluster.mail-type}"  &> dry_run_NAME_DATE.out`

** Main run: **

`snakemake -j 100 --rerun-incomplete --latency-wait 30 --use-conda --cluster-config cluster.yaml --cluster "sbatch -t {cluster.time} --cpus-per-task {cluster.cpus-per-task} --mem {cluster.mem} --mail-type {cluster.mail-type}"  &> main_run_NAME_DATE.out`

In case the pipeline encounters a problem executing the rule 'update_gff_blast', try switching to snakemake v5.19.3 to finish the pipeline, using

`module unload snakemake`

`module load snakemake/5.19.3`


Programs and versions used and tested:


| Program       | Version   |
|---------------|-----------|
| miniconda     | 4.8.3     |
| snakemake     | 6.4.1     |
| snakemake     | 5.19.3    |
| gffread       | 2.2.1     |
| sed           | 4.2.2     |
| python        | 3.7.6     |
| interproscan  | 5.53-87.0 |
| ncbiblastplus | 2.11.0    |
| eggnogmapper  | 2.1.6     |
