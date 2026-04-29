function mutation_ratio_norm = normalize_ratio_by_penta_median(mutation_ratio,penta_idx)
% Normalizes the mutation frequency ratio mutation_ratio (cell array where
% each cell represents a chromosome) by subtracting the pentanucleotide
% medians for the ratio of each pentanucleotide context (the positions of
% all pentanucleotides in the genome is indicated by the cell array
% penta_idx, where each cell is one chromosome. nr_chr represents the
% number of chromosomes in the genome.

nr_chr = length(mutation_ratio);

median_ratio = compute_penta_median_ratio(mutation_ratio,penta_idx);

for chromosome_index = 1:nr_chr
    normalized_ratio = nan(1,length(mutation_ratio{chromosome_index}));
    for pentanucleotide = 1:length(penta_idx{chromosome_index})
        normalized_ratio(penta_idx{chromosome_index}{pentanucleotide}) = mutation_ratio{chromosome_index}(penta_idx{chromosome_index}{pentanucleotide}) - median_ratio(pentanucleotide);
    end
    mutation_ratio_norm{chromosome_index} = normalized_ratio;
end
