function gene_starts = reverse_gene_start_by_strand(annot_genes)
% Computes the gene start coordinates of annot_genes such that the start is
% equal to the stop cooprdinate if the gene is located on the '-'-strand
% (and otherwise the start coordinate is kept).

gene_starts = annot_genes.Start;
for gene_index = 1:length(annot_genes.Start)
    if annot_genes.Strand(gene_index) == '-'
        gene_starts(gene_index) = annot_genes.Stop(gene_index);
    end
end
