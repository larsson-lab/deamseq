function [tri,tri_rv] = find_trinucleotides(ref_genome)
% Computes all trinucleotides of both strands in the genome of ref_genome.

for j = 1:length(ref_genome)
    tri{j} = cell(1,length(ref_genome(j).Sequence));
    tri_rv{j} = cell(1,length(ref_genome(j).Sequence));
    for i = 2:(length(ref_genome(j).Sequence)-1)
        % positions refer to the mutated position
        tri{j}{i} = ref_genome(j).Sequence(i-1:i+1);
        tri_rv{j}{i} = seqrcomplement(tri{j}{i}); 
    end
    tri{j}{1} = '';
    tri{j}{end} = '';
    tri_rv{j}{1} = '';
    tri_rv{j}{end} = '';
end
