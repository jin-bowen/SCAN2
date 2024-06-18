#!/bin/sh
# properties = {"type": "single", "rule": "bwa_mem", "local": false, "input": ["original_fasta/chm13v2.0.fa", "synthetic_reads/chr10/forward_split/00003.fasta", "synthetic_reads/chr10/reverse_split/00003.fasta", "original_fasta/chm13v2.0.fa.amb", "original_fasta/chm13v2.0.fa.ann", "original_fasta/chm13v2.0.fa.bwt", "original_fasta/chm13v2.0.fa.pac", "original_fasta/chm13v2.0.fa.sa"], "output": ["synthetic_reads/chr10/aligned/00003.sam"], "wildcards": {"contig": "chr10", "chunk": "00003"}, "params": {"M": "-M", "K": "-K 100000000"}, "log": ["synthetic_reads/chr10/aligned/00003.log"], "threads": 2, "resources": {"mem_mb": 128000, "mem_mib": 122071, "disk_mb": 22872, "disk_mib": 21813, "tmpdir": "<TBD>"}, "jobid": 71, "cluster": {}}
cd /n/scratch/users/b/boj924/SCAN2-joe && /home/boj924/miniconda3/envs/scan2/bin/python3.10 -m snakemake --snakefile '/home/boj924/scan2_binning/Snakefile' --target-jobs 'bwa_mem:contig=chr10,chunk=00003' --allowed-rules 'bwa_mem' --cores 'all' --attempt 1 --force-use-threads  --resources 'mem_mb=128000' 'mem_mib=122071' 'disk_mb=22872' 'disk_mib=21813' --wait-for-files '/n/scratch/users/b/boj924/SCAN2-joe/.snakemake/tmp.uw_tqq27' 'original_fasta/chm13v2.0.fa' 'synthetic_reads/chr10/forward_split/00003.fasta' 'synthetic_reads/chr10/reverse_split/00003.fasta' 'original_fasta/chm13v2.0.fa.amb' 'original_fasta/chm13v2.0.fa.ann' 'original_fasta/chm13v2.0.fa.bwt' 'original_fasta/chm13v2.0.fa.pac' 'original_fasta/chm13v2.0.fa.sa' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'params' 'mtime' 'code' 'software-env' 'input' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --latency-wait 5 --scheduler 'ilp' --scheduler-solver-path '/home/boj924/miniconda3/envs/scan2/bin' --default-resources 'mem_mb=max(2*input.size_mb, 1000)' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' --mode 2 && touch '/n/scratch/users/b/boj924/SCAN2-joe/.snakemake/tmp.uw_tqq27/71.jobfinished' || (touch '/n/scratch/users/b/boj924/SCAN2-joe/.snakemake/tmp.uw_tqq27/71.jobfailed'; exit 1)

