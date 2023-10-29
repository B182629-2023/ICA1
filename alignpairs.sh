#!/bin/bash

set -e

TCOGENOME="/localdisk/data/BPSM/ICA1/Tcongo_genome/" #Path to reference genome files
FASTQDATA="/localdisk/data/BPSM/ICA1/fastq" #Path to fastq sequence data files
ICA1DIR="/localdisk/data/BPSM/ICA1"


#Creates an index for the reference Tco genome

cd ~/analysis/ref_genome ; cp -r $TCOGENOME/*.gz . 
bowtie2-build -q TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz ref_genome.fa
echo "Reference genome index generated"
cd ~/analysis/seq_reads


#Creates a function which aligns read pairs using bowtie2 and converts to .bam format using samtools

alignreads() {
        sample=$(basename $1 | awk -F"/" '{print $NF;}' | awk -F"[-_]" '{print $2;}') #Basename $1 selects the file name and pipes it to awk to extract the Tco sample ID
        bowtie2 -x ~/analysis/ref_genome/ref_genome.fa -1 ~/analysis/seq_reads/Tco-${sample}_1.fq -2 ~/analysis/seq_reads/Tco-${sample}_2.fq -S ${sample}_alignedreads.sam
	samtools view -bS ${sample}_alignedreads.sam -o ${sample}_alignedreads.bam
	samtools sort ${sample}_alignedreads.bam -o ${sample}_sorted_alignedreads.bam
	samtools index ${sample}_sorted_alignedreads.bam
	bedtools intersect -a $ICA1DIR/TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${sample}_sorted_alignedreads.bam -c > ~/analysis/counts/${sample}_count #-c counts number of overlapping features
}


#Calls the function to all Tco sample fq files

for file in ~/analysis/seq_reads/*1.fq ;
do
        alignreads $file & #This loop takes a while to run, the & symbol runs this in the background
done

wait

echo "Count data for all  Tco samples generated - see ~/analysis/counts/"

#For all Tco counts files, copy file into corresponding Tco sample file in groups directory

for file in ~/analysis/counts/*count; do
        name=$(basename -s _count $file)
	cp $file ~/analysis/groups/*dir/${name}
done

\rm ~/analysis/seq_reads/*.bai
\rm ~/analysis/seq_reads/*.bam


