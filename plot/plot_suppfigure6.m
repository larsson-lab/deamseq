addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat

% Look at pooled sync. / naked ratio:
sync_naked_ratio = mutation_ratio{1}.penta_norm;

% Load ChIP-exo coverage data for Reb1 sites:
Reb1_chip_depth = readtable('../data/ChIP/Reb1_Rap1/SRR346400_depth.tsv','FileType','text','ReadVariableNames',false,'Delimiter','\t');

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


%% Calculate ratio around each predicted TF site
site.ratio = compute_property_around_sites(yeast_genome,sync_naked_ratio,22,site);


%% Re-format the ChIP-exo coverage data
chip_coverage = cell(1,length(yeast_genome));
for chromosome_index = 1:length(yeast_genome)
    chip_coverage{chromosome_index} = zeros(1,length(yeast_genome(chromosome_index).Sequence));
    coverage_current_chromosome = Reb1_chip_depth.Var3(strcmp(string(Reb1_chip_depth.Var1),{yeast_genome(chromosome_index).Header}));
    chip_coverage{chromosome_index}(1:length(yeast_genome(chromosome_index).Sequence)) = coverage_current_chromosome;
end


%% Calculate ChIP-exo coverage around each predicted TF site
site.chip_cov = compute_property_around_sites(yeast_genome,chip_coverage,1000,site);


%% Find correlation with ChIP-exo raw coverage
site_chip_mean = mean(site.chip_cov(990:1011,:),1,'omitnan');
site_ratio_mean = mean(site.ratio(23+3:23+5,:),1,'omitnan');

figure()
plot(site_ratio_mean,log2(site_chip_mean+1),'.')
xlim([-4 2])

[rho,pval] = corr(site_ratio_mean',log2(site_chip_mean+1)','Type','Pearson','Rows','complete');
disp(rho)
disp(pval)
