# E. coli AMR SNP Phylogeny Pipeline

This project performs whole genome analysis of *Escherichia coli* isolates to study antimicrobial resistance (AMR) and evolutionary relationships using SNP phylogeny.

## Workflow

1. Quality control of raw reads
2. Read trimming
3. Genome assembly
4. Assembly quality assessment
5. Genome annotation
6. AMR gene detection
7. Virulence gene detection
8. SNP variant calling
9. Core SNP phylogeny construction

## Tools Used

- FastQC
- MultiQC
- fastp
- QUAST
- Prokka
- Abricate
- Snippy
- IQ-TREE



```pipeline
Raw Reads
   ↓
FastQC / MultiQC
   ↓
fastp (trimming)
   ↓
Genome Assembly
   ↓
QUAST (assembly QC)
   ↓
Prokka (annotation)
   ↓
Abricate (AMR + virulence)
   ↓
Snippy (variant calling)
   ↓
snippy-core (core SNPs)
   ↓
IQ-TREE (phylogeny) 
```
## Pipeline

Run the pipeline:

```bash
bash execute.sh
```
## Outputs
Assembly statistics
Genome annotations
AMR gene profiles
Virulence gene profiles
Core SNP alignment
SNP phylogenetic tree
Phylogenetic Tree

The SNP phylogeny was constructed using IQ-TREE and visualized using iTOL.

Author
Annah Ndono



