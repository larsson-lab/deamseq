function mutation_ratio = compute_mutation_ratio(f,index1,index2)
% Computes the log2 ratio of mutation fractions in f, between 2 samples
% indicated by index1 and index2. nr_chr is the number of chromosomes of
% the genome.
nr_chr = length(f);

mutation_ratio = cell(1,nr_chr);
for chr_idx = 1:nr_chr
    mutation_ratio{chr_idx} = log2(f{chr_idx}(index1,:) ./ f{chr_idx}(index2,:));
    mutation_ratio{chr_idx}(mutation_ratio{chr_idx} == Inf | mutation_ratio{chr_idx} == -Inf) = nan;
end
