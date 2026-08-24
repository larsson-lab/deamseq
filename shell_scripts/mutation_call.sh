# mutation_call.sh <bam_prefix>
# Call mutations using VarScan while retaining all positions, even those lacking substitutions.
# Only relevant columns are retained. Adapt paths below.

reference="/path/to/sacCer3.fa"
varscan="/path/to/VarScan.jar"

samtools mpileup -f $reference -B -q 10 -Q25 -d 1000000 $1.bam | java -jar $varscan mpileup2cns --min-coverage 0 --min-reads2 0 --strand-filter 0 --p-value 1.00 --min-var-freq 0 | cut -f1-5 | tr ':' '\t' | cut -f1-4,6-8 > $1.tsv

