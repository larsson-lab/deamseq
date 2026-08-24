# align.sh <fastq_prefix>
# Aligns, merges and sorts fastq files, with names starting with <fastq_prefix>
# NOTE: assumes paired-end data from two lanes, which will be merged.
# Adapt as needed below, including file suffixes and path to reference genome.

reference="/path/to/sacCer3.fa"

# align
bwa mem -t 15 $reference $1_L001_R1_001.fastq.gz $1_L001_R2_001.fastq.gz | samtools view -bS - > $1_L001.bam &
bwa mem -t 15 $reference $1_L002_R1_001.fastq.gz $1_L002_R2_001.fastq.gz | samtools view -bS - > $1_L002.bam &
wait

# sort
samtools sort -@ 4 -m 4G $1_L001.bam -o $1_L001_sorted.bam &
samtools sort -@ 4 -m 4G $1_L002.bam -o $1_L002_sorted.bam &
wait

# merge
samtools merge $1_sorted.bam $1_L001_sorted.bam $1_L002_sorted.bam

# index
samtools index $1_sorted.bam
