function median_ratio = compute_penta_median_ratio(mutation_ratio,penta_idx)
% Computes the median for mutation_ratio (cell array where each cell
% represents a chromosome) of each pentanucleotide context (penta_idx
% indicates the indices of all pentanucleotides - also cell array where
% each cell is a chromosome). nr_chr indicates the total number of
% chromosomes.

nr_chr = length(mutation_ratio);

median_ratio = nan(1,length(penta_idx{1}));
for pentanucleotide = 1:length(penta_idx{1})
    mutation_ratio_current_pentanucleotide = [];
    for chr_index = 1:nr_chr
        mutation_ratio_current_pentanucleotide = [mutation_ratio_current_pentanucleotide mutation_ratio{chr_index}(penta_idx{chr_index}{pentanucleotide})];
    end

    median_ratio(pentanucleotide) = median(mutation_ratio_current_pentanucleotide,'omitnan');
end
