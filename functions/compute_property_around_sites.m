function ratio_around = compute_property_around_sites(ref_genome,mutation_ratio,plus_minus,site)
% Computes the input property (on the form property{chromosome}(position))
% around sites of ref_genome specified by site.chr, site.start and
% site.stop, modifying the direction based on site.strand. Requires
% site.start <= site.stop.

ratio_around = nan(2*plus_minus+(site.stop(1)-site.start(1)+1),length(site.chr));
for pos_idx = 1:length(site.chr)
    current_chromosome = strcmp({ref_genome.Header},string(site.chr(pos_idx)));

    if site.start(pos_idx) < plus_minus || site.stop(pos_idx) + plus_minus > length(ref_genome(current_chromosome).Sequence) % Move on if gene is too close to edge of chromosome
        continue
    end

    % Reverse order if on reverse complementary strand:
    if site.strand(pos_idx) == '+'
        ratio_around(:,pos_idx) = mutation_ratio{current_chromosome}(1,site.start(pos_idx) - plus_minus:site.stop(pos_idx) + plus_minus);  
    else
        ratio_around(:,pos_idx) = mutation_ratio{current_chromosome}(1,site.stop(pos_idx) + plus_minus:-1:site.start(pos_idx) - plus_minus);
    end
end
