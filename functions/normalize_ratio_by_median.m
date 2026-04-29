function mutation_ratio_norm = normalize_ratio_by_median(mutation_ratio)
% Normalizes the mutation frequency ratio mutation_ratio (cell array where
% each cell represents a chromosome) by subtracting the total median.
% nr_chr represents the number of chromosomes.

nr_chr = length(mutation_ratio);

median_ratio = median(cell2mat(mutation_ratio),'omitnan');
for chr_idx = 1:nr_chr
    mutation_ratio_norm{chr_idx} = mutation_ratio{chr_idx} - median_ratio;
end
