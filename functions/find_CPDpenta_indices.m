function [penta_idx,penta_CPD_compatible] = find_CPDpenta_indices(ref_genome)
% Finds the indices of CPD-compatible pentanucleotides in the whole genome
% of ref_genome, along with the actual pentanucleotide sequences.

trinucleotides_CPD_compatible = {'ACC','ACT','CCA','CCC','CCG','CCT','GCC','GCT','TCA','TCC','TCG','TCT'};
nucleotides = {'A','C','G','T'};

% Construct all possible nucleotides with a central CPD-compatible cytosine:
penta_CPD_compatible = cell(1,4^2*length(trinucleotides_CPD_compatible));
idx = 1;
for tri_idx = 1:length(trinucleotides_CPD_compatible)
    for j = 1:4
        for k = 1:4
            penta_CPD_compatible{idx} = [nucleotides{j} trinucleotides_CPD_compatible{tri_idx} nucleotides{k}];
            idx = idx + 1;
        end
    end
end

% Save index positions for each CPD-pentanucleotide (including reverse complement):
penta_idx = cell(1,length(ref_genome));
for chr_idx = 1:length(ref_genome)
    penta_idx{chr_idx} = cell(1,length(penta_CPD_compatible));
    for pentanucleotide = 1:length(penta_CPD_compatible)
        current_penta = penta_CPD_compatible{pentanucleotide};
        penta_idx{chr_idx}{pentanucleotide} = find_all_pattern_matches(ref_genome(chr_idx).Sequence,current_penta);
        penta_idx{chr_idx}{pentanucleotide} = [penta_idx{chr_idx}{pentanucleotide} find_all_pattern_matches(ref_genome(chr_idx).Sequence,seqrcomplement(current_penta))];
        penta_idx{chr_idx}{pentanucleotide} = penta_idx{chr_idx}{pentanucleotide} + 2;
    end
end
