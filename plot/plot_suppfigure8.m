addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat

% Load nucleosome data:
nucleosomes = readtable('../data/nucleosomes/lift_over/41586_2012_BFnature11142_MOESM263_ESM.reformatted.bed','FileType','text','ReadVariableNames',false,'Delimiter','\t');
nucleosomes = removevars(nucleosomes,"Var3");
nucleosomes = renamevars(nucleosomes,["Var1","Var2","Var5"],["chr","pos","score"]);

% Load motif sequences for transcription factors:
motifs = readtable('../data/motifs/TF_motif_sequences_short.txt','FileType','text','ReadVariableNames',false,'Delimiter',',');
motifs = renamevars(motifs,["Var1","Var2","Var3","Var4"],["Name","Seq","Sym_seq","Sym_seq_rcomp"]);


%% Find positions for Reb1 based on motif sequence
TF_name = 'Reb1';
TF_motif_index = find(strcmp(motifs.Name,TF_name));
motif_sequence = motifs.Seq{TF_motif_index};

% Find motif positions:
site = find_seq_motif_sites(yeast_genome,motif_sequence,motifs.Sym_seq{TF_motif_index});

fprintf(['Number of motif sites for ' TF_name ': %d\n'],length(site.chr));


%% Remove all sites that are not included in the whitelisted regions
outside_whitelisted = false(1,length(site.chr));
for site_index = 1:length(site.chr)
    current_chr = strcmp(site.chr(site_index),{yeast_genome.Header});
    outside_whitelisted(site_index) = any(idx_cov_ok{current_chr}(site.start(site_index):site.stop(site_index)) == 0);
end
site.chr = site.chr(~outside_whitelisted);
site.start = site.start(~outside_whitelisted);
site.stop = site.stop(~outside_whitelisted);
site.strand = site.strand(~outside_whitelisted);
site.motif = site.motif(~outside_whitelisted);


%% Find which nucleosomes have Reb1 binding on either side
nearby_reb1 = false(1,length(nucleosomes.chr));
for chr_index = 1:length(yeast_genome)
    nucleosome_index = find(strcmp(nucleosomes.chr,string(yeast_genome(chr_index).Header)));

    current_site_start = site.start(strcmp(site.chr,string(yeast_genome(chr_index).Header)));
    current_site_stop = site.stop(strcmp(site.chr,string(yeast_genome(chr_index).Header)));
    for site_index = 1:length(current_site_start)
        idx = nucleosomes.pos(nucleosome_index) - 200 < current_site_start(site_index) & nucleosomes.pos(nucleosome_index) + 200 > current_site_stop(site_index);

        % If current Reb1 site is bound, set nearby_reb1=true for all nucleosomes within +-200 bp:
        if max(abs(mutation_ratio{4}.penta_norm{chr_index}(current_site_start(site_index):current_site_stop(site_index)))) >= 2
            nearby_reb1(nucleosome_index(idx)) = true;
        end
    end
end


%% Mask the ratio at Reb1-sites
mutation_ratio_masked = mutation_ratio{4}.penta_norm;
mutation_ratio_aa_masked = mutation_ratio{5}.penta_norm;
for chr_index = 1:length(yeast_genome)
    current_site_start = site.start(strcmp(site.chr,string(yeast_genome(chr_index).Header)));
    current_site_stop = site.stop(strcmp(site.chr,string(yeast_genome(chr_index).Header)));
    sites = [];
    for site_index = 1:length(current_site_start)
        sites = [sites current_site_start(site_index)-5:current_site_stop(site_index)+5];
    end

    mutation_ratio_masked{chr_index}(sites) = nan;
    mutation_ratio_aa_masked{chr_index}(sites) = nan;
end


%% Calculate ratio relative nucleosomes (sorted after nucleosome occupancy)
nucleosomes.start = nucleosomes.pos;
nucleosomes.stop = nucleosomes.pos;
nucleosomes.strand = strings(1,length(nucleosomes.chr))';
nucleosomes.strand(:) = '+';

[~,sorting_index] = sort(nucleosomes.score,'descend','MissingPlacement','last');

nucleosomes.chr = nucleosomes.chr(sorting_index);
nucleosomes.start = nucleosomes.start(sorting_index);
nucleosomes.stop = nucleosomes.stop(sorting_index);
nucleosomes.score = nucleosomes.score(sorting_index);
nucleosomes.strand = nucleosomes.strand(sorting_index);


ratio_nucleosome = compute_property_around_sites(yeast_genome,mutation_ratio_masked,1000,nucleosomes);
ratio_nucleosome_rap = compute_property_around_sites(yeast_genome,mutation_ratio_aa_masked,1000,nucleosomes);

nearby_reb1 = nearby_reb1(sorting_index);


%% Plot mean ratio relative nucleosome positions for nucleosomes with Reb1 sites nearby
xpos = -1000:1000;
figure()

subplot(2,2,1)
plot(xpos,mean(ratio_nucleosome(:,nearby_reb1),2,'omitnan'),LineWidth=1);
xlim([-200 200])
ylim([-0.2 0.15])
title(['Reb1 sites, -rap, n=' num2str(length(find(nearby_reb1)))])

subplot(2,2,3)
plot(xpos,mean(ratio_nucleosome_rap(:,nearby_reb1),2,'omitnan'),LineWidth=1);
xlim([-200 200])
ylim([-0.2 0.15])
title(['Reb1 sites, +rap, n=' num2str(length(find(nearby_reb1)))])

subplot(2,2,2)
plot(xpos,mean(ratio_nucleosome(:,~nearby_reb1),2,'omitnan'),LineWidth=1);
xlim([-200 200])
ylim([-0.2 0.15])
title(['No Reb1 sites, -rap, n=' num2str(length(find(~nearby_reb1)))])

subplot(2,2,4)
plot(xpos,mean(ratio_nucleosome_rap(:,~nearby_reb1),2,'omitnan'),LineWidth=1);
xlim([-200 200])
ylim([-0.2 0.15])
title(['No Reb1 sites, +rap, n=' num2str(length(find(~nearby_reb1)))])

