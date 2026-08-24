# align.sh <fastq_prefix>
# Aligns, merges and sorts fastq files, with names starting with <fastq_prefix>
# NOTE: assumes single-end data from one lane.
# Adapt as needed below, including file suffixes and path to reference genome.

reference="/path/to/sacCer3.fa"

# align
bwa mem -t 15 $reference $1.fastq.gz | samtools view -bS - > $1.bam &
wait

# sort
samtools sort -@ 4 -m 4G $1.bam -o $1_sorted.bam &
wait

# index
samtools index $1_sorted.bam

# Calculate per base coverage
samtools depth -a $1_sorted.bam > $1_depth.tsv
