<img src="https://github.com/user-attachments/assets/c1778e95-639c-494f-8550-552cde5c5a8e" alt="pipeline" width="300" height="337.5">

## 1. Dependent Software

- Java
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [SignalP6](https://services.healthtech.dtu.dk/services/SignalP-6.0/)
- [Predisi](http://predisi.de/)
- [Phobius](https://phobius.sbc.su.se/)

## 2. What to input

- If the primary transcript's protein sequence FASTA file is provided, add the protein file.
- If the primary transcript's protein sequence FASTA file is not provided, add the genome file and the annotation file.

gzipped files are supported.

## 3. What to output

- Signal peptides prediction results from three softwares
- Final phytocytokine candidates

## 4. Usage

```shell
snakemake \
	--snakefile SecProScan.smk \
	--use-conda \
	--use-singularity \
	--rerun-incomplete \
	--nolock
```
