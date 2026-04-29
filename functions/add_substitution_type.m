function mut = add_substitution_type(ref_genome,mut)
% Adds the type of substitution of each mutation as a new field (.subst)
% in the struct mut (requires a field mut.var) for the nucleotide in the
% mutated data. ref_genome is a .fa file with the reference genome.

n_files = size(mut.n{1},1);
for chr_index = 1:length(ref_genome)
    mut.subst{chr_index} = zeros(n_files, length(ref_genome(chr_index).Sequence), 'int8'); % substitution type: 0 = none, 1-6 = C>A, C>G, C>T, T>A, T>C, T>G, 7=indel
    for sample_index = 1:n_files
        mut.subst{chr_index}(sample_index,(ref_genome(chr_index).Sequence == 'C' & mut.var{chr_index}(sample_index,:) == 'A') | (ref_genome(chr_index).Sequence == 'G' & mut.var{chr_index}(sample_index,:) == 'T')) = 1;
        mut.subst{chr_index}(sample_index,(ref_genome(chr_index).Sequence == 'C' & mut.var{chr_index}(sample_index,:) == 'G') | (ref_genome(chr_index).Sequence == 'G' & mut.var{chr_index}(sample_index,:) == 'C')) = 2;
        mut.subst{chr_index}(sample_index,(ref_genome(chr_index).Sequence == 'C' & mut.var{chr_index}(sample_index,:) == 'T') | (ref_genome(chr_index).Sequence == 'G' & mut.var{chr_index}(sample_index,:) == 'A')) = 3;
        mut.subst{chr_index}(sample_index,(ref_genome(chr_index).Sequence == 'T' & mut.var{chr_index}(sample_index,:) == 'A') | (ref_genome(chr_index).Sequence == 'A' & mut.var{chr_index}(sample_index,:) == 'T')) = 4;
        mut.subst{chr_index}(sample_index,(ref_genome(chr_index).Sequence == 'T' & mut.var{chr_index}(sample_index,:) == 'C') | (ref_genome(chr_index).Sequence == 'A' & mut.var{chr_index}(sample_index,:) == 'G')) = 5;
        mut.subst{chr_index}(sample_index,(ref_genome(chr_index).Sequence == 'T' & mut.var{chr_index}(sample_index,:) == 'G') | (ref_genome(chr_index).Sequence == 'A' & mut.var{chr_index}(sample_index,:) == 'C')) = 6;
        mut.subst{chr_index}(sample_index,(mut.var{chr_index}(sample_index,:) == '+' | mut.var{chr_index}(sample_index,:) == '-')) = 7;
    end
end
