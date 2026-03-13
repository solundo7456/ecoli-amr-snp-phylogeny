#!/bin/bash

###############################
#first first qc
###############################

base='/home/ndono_school/Downloads/bioinformatics/AMR/clin'

cd ${base}/raw || exit

mkdir -p ${base}/qc 
echo 'Doing fastqc'
fastqc *.gz -o ${base}/qc 

mkdir -p ${base}/multiqc
echo'Running multiqc'
multiqc ${base}/qc  -o ${base}/multiqc 


###############################
# trimming 
###############################
mkdir -p ${base}/trimmed
cd ${base}/raw || exit

for f in *fastq.gz; do
	
	basen=$(basename "$f" .fastq.gz)
	echo "Processing $basen"
	fastp -i "$f" \
	      -o "${base}/trimmed/${basen}_trimmed.fastq.gz" \
	      -w 8 \
	      --trim_front1 15
done
###############################
# fastqc trimmed
###############################
mkdir -p ${base}/trimmed_qc
mkdir -p ${base}/trimmed_multiqc

cd ${base}/trimmed || exit
 
echo 'Doing fastqc'
fastqc *.gz -o ${base}/trimmed_qc 


echo 'Running multiqc'
multiqc ${base}/trimmed_qc  -o ${base}/trimmed_multiqc 

#openning html
###############################
cd ${base}/trimmed_qc || exit
for i in *.html; do 
	xdg-open ${i}
done

cd ${base}/trimmed_multiqc || exit
for i in *.html; do
	xdg-open ${i}
done


###############################
#assembly
################################


# Trimmed assemblies
################################

mkdir -p ${base}/quast_trimmed

for i in ${base}/trim_assemble/*/contigs.fasta; do

    sample=$(basename $(dirname "$i"))
    echo "Running QUAST for trimmed sample $sample"

    quast.py "$i" \
    -o ${base}/quast_trimmed/${sample} \
    -t 8

done



# Raw assemblies
################################

mkdir -p ${base}/quast_raw

for i in ${base}/raw_assemble/*/contigs.fasta; do

    sample=$(basename $(dirname "$i"))
    echo "Running QUAST for raw sample $sample"

    quast.py "$i" \
    -o ${base}/quast_raw/${sample} \
    -t 8

done


# Compare raw and trimmed assemblies
################################

quast.py \
trim_assemble/*/contigs.fasta \
raw_assemble/*/contigs.fasta \
-o quast_compare \
-t 8



########################################
#assembly qc
########################################

threads=8

# QUAST for trimmed assemblies
########################################

mkdir -p ${base}/quast_trimmed

for i in ${base}/trim_assemble/*/contigs.fasta; do

    sample=$(basename "$(dirname "$i")")
    echo "Running QUAST for trimmed sample $sample"

    quast.py "$i" \
        -o ${base}/quast_trimmed/${sample} \
        -t $threads

done


# QUAST for raw assemblies
########################################

mkdir -p ${base}/quast_raw

for i in ${base}/raw_assemble/*/contigs.fasta; do

    sample=$(basename "$(dirname "$i")")
    echo "Running QUAST for raw sample $sample"

    quast.py "$i" \
        -o ${base}/quast_raw/${sample} \
        -t $threads

done



# QUAST comparison
########################################

mkdir -p ${base}/quast_compare

echo "Running QUAST comparison"

quast.py \
${base}/trim_assemble/*/contigs.fasta \
${base}/raw_assemble/*/contigs.fasta \
-o ${base}/quast_compare \
-t $threads

echo "QUAST analysis finished"



################################
# genome annotation
################################

mkdir -p ${base}/annotation

for i in ${base}/trim_assemble/*/contigs.fasta;do
	basen=$(basename (dirname "$i"))
	echo "Running Genome annotaion on $basen"
	prokka "$i" --outdir ${base}/annotation/${basen} \
	--prefix "${basen}" --genus Escherichia \
	--species coli --strain O157 \
	--cpus 8 --usegenus --rfam
done 


################################
# AMR using abricate 
################################

mkdir -p ${base}/abricate/amr_abricate

for i in ${base}/trim_assemble/*/contigs.fasta;do
	basen=$(basename $(dirname "$i"))
	echo " Running abricate amr scan on $basen"
	for db in resfinder card; do
		abricate "$i" --db "${db}" > ${base}/abricate/amr_abricate/${basen}_${db}.txt
	done			
done

# AMR using abricate summary
################################
echo "Generating combined AMR summary..."
abricate --summary ${base}/abricate/amr_abricate/*.txt > ${base}/abricate/amr_abricate/amr_abricate_summary.txt

echo "Done! Summary file is at: ${base}/amr_abricate/amr_abricate_summary.txt"


# virulence using abricate 
################################
mkdir -p ${base}/abricate/virulence_abricate


for i in ${base}/trim_assemble/*/contigs.fasta;do
	basen=$(basename $(dirname "$i"))
	echo " Running abricate virulence scan on $basen"
	abricate --db vfdb "$i" > ${base}/abricate/virulence_abricate/${basen}.txt
done


# AMR using abricate summary
################################
echo "Generating combined virulence summary..."
abricate --summary ${base}/abricate/virulence_abricate/*.txt > ${base}/abricate/virulence_abricate/virulence_abricate_summary.txt

echo "Done! Summary file is at: ${base}/abricate/virulence_abricate/amr_abricate_summary.txt"	


################################
# snippy- variant calling
################################
mkdir -p ${base}/snps
ref=/home/ndono_school/Downloads/bioinformatics/AMR/clin/reference/ecoli_ref.fasta
for i in ${base}/trim_assemble/*/contigs.fasta;do
	basen=$(basename $(dirname "$i"))
	echo "Doing variant calling for $basen"
	snippy --ctgs "$i" --ref "$ref" --outdir ${base}/snps/${basen} --cpus 4
	
done


# variant consolidation
################################
cd ${base}/snps

snippy-core --ref /home/ndono_school/Downloads/bioinformatics/AMR/clin/reference/ecoli_ref.fasta */

################################
# SNPS phylogenetic tree construction
################################

iqtree -s core.aln -m GTR+G -bb 1000 -nt AUTO

