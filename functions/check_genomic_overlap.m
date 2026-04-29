function overlap_index = check_genomic_overlap(chromosomes,region_set_1,region_set_2)
% Check which sites in the struct region_set_1 (requires fields start,
% stop and chr) that overlaps with those in struct region_set_2 (also
% requires fields start, stop and chr). chromosomes represents the
% annotation used for the chromosomes names.

overlap_index = zeros(1,length((region_set_1.chr)));
for site_index = 1:length(region_set_1.chr)
    current_chromosome = strcmp(region_set_1.chr(site_index),string(chromosomes));
    current_site_start = min(region_set_1.start(site_index),region_set_1.stop(site_index));
    current_site_end = max(region_set_1.start(site_index),region_set_1.stop(site_index));
    current_site_coords = current_site_start:current_site_end;

    footprints_start_current_chromosome = region_set_2.start(strcmp(region_set_2.chr,chromosomes{current_chromosome}));
    footprints_end_current_chromosome = region_set_2.stop(strcmp(region_set_2.chr,chromosomes{current_chromosome}));

    for footprint_index = 1:length(footprints_start_current_chromosome)
        footprint_coord = footprints_start_current_chromosome(footprint_index):footprints_end_current_chromosome(footprint_index);
        if any(ismember(current_site_coords,footprint_coord))
            overlap_index(site_index) = 1;
        end
    end
end
