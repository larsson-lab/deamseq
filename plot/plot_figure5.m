addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat

% Load binding data for sites:
Reb1_chip = readtable('../data/ChIP/Reb1_Rap1/lift_over/sacCer3.Reb1_peaks.bed','FileType','text','ReadVariableNames',false);

% Load motif sequences for transcription factors:
motifs = readtable('../data/motifs/TF_motif_sequences_short.txt','FileType','text','ReadVariableNames',false,'Delimiter',',');
motifs = renamevars(motifs,["Var1","Var2","Var3","Var4"],["Name","Seq","Sym_seq","Sym_seq_rcomp"]);

% Load trans. reg. code:
pred_sites = readtable('../data/transRegCode/transRegCode.txt','FileType','text','ReadVariableNames',false,'Delimiter','tab');
pred_sites = removevars(pred_sites,["Var1","Var7","Var8"]);
pred_sites = renamevars(pred_sites,["Var2","Var3","Var4","Var5","Var6"],["chr","start","stop","name","score"]);
pred_sites.start = pred_sites.start + 1;


%% Find the positions of strongest difference and look for meme motif
% Extract maximum abs. signal of each bin, and exact genomic position of that signal:
bin_size = 10;
binned_signal = cell(1,length(yeast_genome));
binned_signal_pos = cell(1,length(yeast_genome)); % Pos of strongest signal in bin
for chromosome_index = 1:length(yeast_genome)
    ratio_diff = abs(mutation_ratio{4}.penta_norm{chromosome_index} - mutation_ratio{5}.penta_norm{chromosome_index});

    edited_ratio_diff = [ratio_diff nan(1,bin_size-mod(length(mutation_ratio{4}.penta_norm{chromosome_index}),bin_size))];
    nr_bins = length(edited_ratio_diff) / bin_size;
    ratio_matrix = reshape(edited_ratio_diff,[bin_size,nr_bins]);

    for bin_index = 1:nr_bins
        binned_signal{chromosome_index}(bin_index) = max(ratio_matrix(:,bin_index));
        index_max = ratio_matrix(:,bin_index) == max(ratio_matrix(:,bin_index));

        index = nan;
        if isscalar(find(index_max))
            index = bin_size * (bin_index - 1) + find(index_max) - 1;
        end

        binned_signal_pos{chromosome_index}(bin_index) = index;
    end
end


% Iterate through the data and extract the strongest signals:
strongest_pos = [];
strongest_chr = [];
strongest_signal = [];
for chromosome_index = 1:length(yeast_genome)
    strongest_idx = find(binned_signal{chromosome_index} >= 2);

    strongest_signal = [strongest_signal binned_signal{chromosome_index}(strongest_idx)];
    strongest_pos = [strongest_pos binned_signal_pos{chromosome_index}(strongest_idx)];
    strongest_chr = [strongest_chr chromosome_index * ones(1,length(strongest_idx))];
end

% Only keep the 100 strongest signals
[~,sorting_index] = sort(strongest_signal,'descend');
strongest_pos = strongest_pos(sorting_index(1:100));
strongest_chr = strongest_chr(sorting_index(1:100));

% Find sequence contexts around the 100 strongest signals, to run through MEME motif discovery:
seq_context = strings([]);
for index = 1:length(strongest_pos)
    seq_context(index) = yeast_genome(strongest_chr(index)).Sequence(strongest_pos(index)-7:strongest_pos(index)+7);
end

writematrix(seq_context,'fig5_seqcontext.txt')


%% Remove sites in trans. reg. code dataset that are partly/completely outside whitelisted regions
outside_whitelisted = false(1,length(pred_sites.chr));
for site_index = 1:length(pred_sites.chr)
    current_chr = strcmp(pred_sites.chr(site_index),{yeast_genome.Header});
    current_site_coords = pred_sites.start(site_index):pred_sites.stop(site_index);
    
    idx_cov_ok_current_chromosome = idx_cov_ok{current_chr};
    outside_whitelisted(site_index) = any(idx_cov_ok_current_chromosome(current_site_coords) == 0);
end

pred_sites = pred_sites(~outside_whitelisted,:);


%% Calculate DNAse I footprint overlap for each site
pred_sites.footprint = check_genomic_overlap({yeast_genome.Header},pred_sites,dnase_footprints)';


%% Calculate mean abs ratios (Reb1aa and Reb1aa+rap) for each site
pred_sites.ratio_aa = nan(length(pred_sites.chr),1);
pred_sites.ratio_aa_rap = nan(length(pred_sites.chr),1);
for site_index = 1:length(pred_sites.chr)
    current_chr = strcmp(pred_sites.chr(site_index),{yeast_genome.Header});
    current_site_coords  = pred_sites.start(site_index):pred_sites.stop(site_index);

    pred_sites.ratio_aa(site_index) = mean(abs(mutation_ratio{4}.penta_norm{current_chr}(current_site_coords)),'omitnan');
    pred_sites.ratio_aa_rap(site_index) = mean(abs(mutation_ratio{5}.penta_norm{current_chr}(current_site_coords)),'omitnan');
end


%% Wilcoxon signrank test (with Bonferroni correction) of all trans. reg. code. TFs to compare the untreated and treated
TF_array = unique(pred_sites.name);

q = nan(1,length(TF_array));
ratio_diff = nan(1,length(TF_array));
for TF_index = 1:length(TF_array)
    TF_sites = find(strcmp(pred_sites.name,TF_array(TF_index)));
    dnase_overlap = pred_sites.footprint(TF_sites);
    ratio_aa = pred_sites.ratio_aa(TF_sites);
    ratio_aa_rap = pred_sites.ratio_aa_rap(TF_sites);

    TF_score_ratios{1} = ratio_aa;
    TF_score_ratios{2} = ratio_aa_rap;

    % Wilcoxon rank sum test with Bonferroni correction of p-values:
    if ~isempty(TF_score_ratios{1}) && ~isempty(TF_score_ratios{2})
        rank_sum = signrank(TF_score_ratios{1},TF_score_ratios{2});
        q(TF_index) = rank_sum * length(TF_array); % Bonferroni correction
        ratio_diff(TF_index) = mean(TF_score_ratios{1},'omitnan') - mean(TF_score_ratios{2},'omitnan');
    end
end

q(q > 1) = 1;


%% Bonferroni-correction of p-values and plot results
figure()
subplot(3,4,[1 2;5 6])
hold on

h = line([0 120],[0 0],'color','k');
plot(-log10(q),ratio_diff,'o',MarkerSize=8,MarkerFaceColor=[0.85 0.85 0.85],MarkerEdgeColor='none');
xline(-log10(0.01),LineStyle="--");

[sorted_q,index_sorted] = sort(q,'ascend');
sorted_ratio_diff = ratio_diff(index_sorted);
sorted_TF_array = TF_array(index_sorted);
plot(-log10(sorted_q(1:length(find(q <= 0.01)))),sorted_ratio_diff(1:length(find(q <= 0.01))),'o',MarkerSize=8);
text(-log10(sorted_q(1:length(find(q <= 0.01))))+1,sorted_ratio_diff(1:length(find(q <= 0.01))),sorted_TF_array(1:length(find(q <= 0.01))));

xlabel('-log10(p)')
ylabel('Ratio diff')
xlim([0 120])
uistack(h,'bottom');


%% Save supplementary data 2
P_value = q';
Ratio_difference = -ratio_diff';
T = table(TF_array,P_value,Ratio_difference);
writetable(T,'supplementary_data_2.xls','FileType','spreadsheet')


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
site.ratio_aa = compute_property_around_sites(yeast_genome,mutation_ratio{4}.penta_norm,22,site);
site.ratio_aa_rap = compute_property_around_sites(yeast_genome,mutation_ratio{5}.penta_norm,22,site);


%% Investigate if motif sites overlap with DNase I footprints
site.footprint = check_genomic_overlap({yeast_genome.Header},site,dnase_footprints);


%% Check overlap with ChIP-peaks
if strcmp(TF_name,'Reb1')
    chip.chr = Reb1_chip.Var1;
    chip.peak = Reb1_chip.Var2;
    chip.width = Reb1_chip.Var4;
    chip.start = chip.peak - ceil(chip.width / 2);
    chip.stop = chip.peak + ceil(chip.width / 2);

    site.chip = check_genomic_overlap({yeast_genome.Header},site,chip);
end


%% Sort sites based on ratio
site.site_ratio = mean(abs(site.ratio_aa(23:23+length(char(site.motif(1)))-1,:)),1,'omitnan');
[~,sorting_index] = sort(abs(site.site_ratio),'descend','MissingPlacement','last');


fields = fieldnames(site);
for i = 1:length(fieldnames(site))
    site.(fields{i}) = site.(fields{i})(:,sorting_index);
end


%% Plot heatmap
motif_length = length(char(site.motif(1)));
color_map = readmatrix('../colormaps/diverging_colormap_redblue_new2','FileType','text');

figure()
subplot(1,8,[1 3])
imAlpha=ones(size(site.ratio_aa'));
imAlpha(isnan(site.ratio_aa')) = 0;
image(site.ratio_aa',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.8 0.8 0.8]);
colormap(flip(color_map));
clim([-3 3]);
colorbar(Location='westoutside');
sgtitle('Log_2(f_{sync}/f_{naked}) for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:25:51;
ax.XTickLabel = -25:25:25;
xlim([23-22-0.5 23+motif_length-1+22+0.5])

subplot(1,8,[4 6])
imAlpha=ones(size(site.ratio_aa_rap'));
imAlpha(isnan(site.ratio_aa_rap')) = 0;
image(site.ratio_aa_rap',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.8 0.8 0.8]);
colormap(flip(color_map));
clim([-3 3]);
colorbar(Location='westoutside');
sgtitle('Log_2(f_{sync}/f_{naked}) for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:25:51;
ax.XTickLabel = -25:25:25;
xlim([23-22-0.5 23+motif_length-1+22+0.5])

subplot(1,8,7)
barh(flip(site.footprint'),Facecolor=[1 0 0]);
xlabel('Footprints')

subplot(1,8,8)
barh(flip(site.chip'),Facecolor=[0 1 0]);
xlabel('ChIP')

