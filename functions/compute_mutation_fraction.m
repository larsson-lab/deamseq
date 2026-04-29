function f = compute_mutation_fraction(mut)
% Calculates the fraction between number of mutated reads and total number
% of reads specified by the struct (requires fields .n for mutated reads
% and .cov for coverage), with one pseudocount.

nr_chr = length(mut.n);
for chr_index = 1:nr_chr
    f{chr_index} = double(mut.n{chr_index} + 1)./double(mut.cov{chr_index} + 1); % add one pseudocount
end