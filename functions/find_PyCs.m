function idx_c_dipy = find_PyCs(ref_genome)
% Make cell arrays describing the trinucleotide at each position in the
% reference genome ref_genome.

[tri,tri_rv] = find_trinucleotides(ref_genome);

c_dipy = {'CCC' 'CCG' 'CCA' 'CCT' 'TCC' 'TCA' 'TCG' 'ACC' 'GCC' 'ACT' 'GCT' 'TCT'};
for chr_index = 1:length(ref_genome)
    idx_c_dipy{chr_index} = ismember(tri{chr_index}, c_dipy) | ismember(tri_rv{chr_index}, c_dipy);
end
