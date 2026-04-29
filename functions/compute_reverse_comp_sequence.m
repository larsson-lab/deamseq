function motif_rcomp_sequence = compute_reverse_comp_sequence(motif_sequence)
% Computes the reverse complementary sequence of the string motif_sequence,
% with characters compatible with regexp (. meaning any character, []
% indicating OR of the characters within it.

k = strfind(motif_sequence,'.');
l = strfind(motif_sequence,'[');
m = strfind(motif_sequence,']');
index_letters = setdiff(1:length(motif_sequence),[k l m]);
indices = [index_letters [k l m]];
characters = [seqcomplement(motif_sequence(index_letters)) motif_sequence([k m l])]; % Flip position of [ and ]
new_indices = length(motif_sequence) - (indices-1);
[~,idx] = sort(new_indices);
motif_rcomp_sequence = characters(idx);
