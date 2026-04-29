function site = find_seq_motif_sites(ref_genome,motif_sequence,sym_motif_sequence)
% Finds all sites in the reference genome ref_genome which matches the
% sequence pattern given by the regexp-compatible string motif_sequence
% (with symbolic sequence indicated by sym_motif_sequence). Returns a
% struct site with fields .chr, .start. stop indicating the coordinates of
% a match, .strand indicating the strannd location of the match ('+' or
% '-') and .motif indicating the exact motif sequence of each site (since
% the input motif_sequence can contain degenerate positions).

motif_rcomp_sequence = compute_reverse_comp_sequence(motif_sequence);

site.chr = string([]);
site.start = [];
site.stop = [];
site.strand = string([]);
site.motif = string([]);
site_forward = [];
for current_chromosome = 1:length(ref_genome)
    forward_motif_start_positions = find_all_pattern_matches(ref_genome(current_chromosome).Sequence,motif_sequence);
    reverse_motif_start_positions = find_all_pattern_matches(ref_genome(current_chromosome).Sequence,motif_rcomp_sequence);

    % Remove for palindromes:
    reverse_motif_start_positions = reverse_motif_start_positions(~ismember(reverse_motif_start_positions,forward_motif_start_positions));

    % Save:
    site.start = [site.start forward_motif_start_positions reverse_motif_start_positions+(length(sym_motif_sequence)-1)];
    site.stop = [site.stop forward_motif_start_positions+(length(sym_motif_sequence)-1) reverse_motif_start_positions];
    site_forward = [site_forward ones(1,length(forward_motif_start_positions)) zeros(1,length(reverse_motif_start_positions))];

    current_motif_start_positions = [forward_motif_start_positions reverse_motif_start_positions];
    current_motif_chromosome = strings(1,length(current_motif_start_positions));
    current_motif_chromosome(:) = ref_genome(current_chromosome).Header;
    site.chr = [site.chr current_motif_chromosome];

    current_individual_motifs = strings(1,length(current_motif_start_positions));
    for i = 1:length(current_motif_start_positions)
        try
            current_individual_motifs(i) = ref_genome(current_chromosome).Sequence(current_motif_start_positions(i):current_motif_start_positions(i)+(length(sym_motif_sequence)-1));
        end
    end
    site.motif = [site.motif current_individual_motifs];

    for site_index = 1:length(site.start)
        min_coord = min(site.start(site_index),site.stop(site_index));
        max_coord = max(site.start(site_index),site.stop(site_index));
        site.start(site_index) = min_coord;
        site.stop(site_index) = max_coord;
    end
end

% TODO: Gör detta direkt:
site.strand(site_forward == 1) = '+';
site.strand(site_forward == 0) = '-';
