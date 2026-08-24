# dedup.sh <bam_prefix>
# Perform UMI deduplication using UMItools. UMI is assumed to be last in the read name, separated by ':'.
# <bam_prefix> is the bam file name without '.bam'.

umi_tools dedup -I $1.bam --paired --umi-separator=: --chimeric-pairs=discard --unpaired-reads=discard -S $1_dedup.bam
