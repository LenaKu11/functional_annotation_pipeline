Dry mode:

snakemake -npr -j 100 --use-conda --cluster-config cluster.yaml --cluster "sbatch -t {cluster.time} --cpus-per-task {cluster.cpus-per-task} --mem {cluster.mem} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}"  &> dry_run_Cmulti_DATE.out

Main run:

snakemake -j 100 --rerun-incomplete --latency-wait 30 --use-conda --cluster-config cluster.yaml --cluster "sbatch -t {cluster.time} --cpus-per-task {cluster.cpus-per-task} --mem {cluster.mem} --mail-type {cluster.mail-type} --mail-user {cluster.mail-user}"  &> main_run_Cmulti_DATE.out
# functional_annotation_pipeline
