addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat

% Look at pooled sync. / naked ratio:
cell_naked_ratio = mutation_ratio{1}.penta_norm;

% Load expression data:
expression = readtable('../data/expression/FPKM_NR_CR.txt','FileType','text','ReadVariableNames',true,'Delimiter','\t');

% Load binding data for sites:
Reb1_chip = readtable('../data/ChIP/Reb1_Rap1/lift_over/sacCer3.Reb1_peaks.bed','FileType','text','ReadVariableNames',false);
Reb1_chip_depth = readtable('../data/ChIP/Reb1_Rap1/SRR346400_depth.tsv','FileType','text','ReadVariableNames',false,'Delimiter','\t');

% Load motif sequences for transcription factors:
motifs = readtable('../data/motifs/TF_motif_sequences_short.txt','FileType','text','ReadVariableNames',false,'Delimiter',',');
motifs = renamevars(motifs,["Var1","Var2","Var3","Var4"],["Name","Seq","Sym_seq","Sym_seq_rcomp"]);

% Load trans. reg. code:
pred_sites = readtable('../data/transRegCode/transRegCode.txt','FileType','text','ReadVariableNames',false,'Delimiter','tab');
pred_sites = removevars(pred_sites,["Var1","Var7","Var8"]);
pred_sites = renamevars(pred_sites,["Var2","Var3","Var4","Var5","Var6"],["chr","start","stop","name","score"]);
pred_sites.start = pred_sites.start + 1;


%% Remove sites in trans. reg. code dataset that are partly or completely outside whitelisted regions
outside_whitelisted = false(1,length(pred_sites.chr));
for site_index = 1:length(pred_sites.chr)
    current_chr = strcmp(pred_sites.chr(site_index),{yeast_genome.Header});
    current_site_coords = pred_sites.start(site_index):pred_sites.stop(site_index);
    
    idx_cov_ok_current_chromosome = idx_cov_ok{current_chr};
    outside_whitelisted(site_index) = any(idx_cov_ok_current_chromosome(current_site_coords) == 0);
end
pred_sites = pred_sites(~outside_whitelisted,:);


%% Calculate mean abs ratio as well as DNAse I footprint overlap for each site:
pred_sites.dnase_overlap = check_genomic_overlap({yeast_genome.Header},pred_sites,dnase_footprints)';

pred_sites.mean_abs_ratio = nan(length(pred_sites.chr),1);
for site_index = 1:length(pred_sites.chr)
    current_chr = strcmp(pred_sites.chr(site_index),{yeast_genome.Header});
    pred_sites.mean_abs_ratio(site_index) = mean(abs(cell_naked_ratio{current_chr}(pred_sites.start(site_index):pred_sites.stop(site_index))),'omitnan');
end


%% Wilcoxon rank sum test (with Bonferroni correction) of all trans. reg. code. TFs to compare DNase supported and non-supported sites
TF_array = unique(pred_sites.name);

q = nan(1,length(TF_array));
ratio_diff = nan(1,length(TF_array));
for TF_index = 1:length(TF_array)
    TF_sites = find(strcmp(pred_sites.name,TF_array(TF_index)));
    mean_abs_ratio = pred_sites.mean_abs_ratio(TF_sites);
    dnase_overlap = pred_sites.dnase_overlap(TF_sites);

    % Wilcoxon rank sum test with Bonferroni correction of p-values:
    TF_score_ratios{1} = mean_abs_ratio(dnase_overlap == 0);
    TF_score_ratios{2} = mean_abs_ratio(dnase_overlap == 1);
    if ~isempty(TF_score_ratios{1}) && ~isempty(TF_score_ratios{2})
        rank_sum = ranksum(TF_score_ratios{1},TF_score_ratios{2});
        q(TF_index) = rank_sum * length(TF_array); % Bonferroni correction
        ratio_diff(TF_index) = mean(TF_score_ratios{1},'omitnan') - mean(TF_score_ratios{2},'omitnan');
    end
end
q(q > 1) = 1;


%% Plot results
figure()
subplot(3,4,[1 2;5 6])
hold on
plot(-log10(q),-ratio_diff,'o',MarkerSize=8,MarkerFaceColor=[0.85 0.85 0.85],MarkerEdgeColor='none');
xline(-log10(0.01),LineStyle="--");

[sorted_q,index_sorted] = sort(q,'ascend');
sorted_ratio_diff = ratio_diff(index_sorted);
sorted_TF_array = TF_array(index_sorted);
plot(-log10(sorted_q(1:6)),-sorted_ratio_diff(1:6),'o',MarkerSize=8);
text(-log10(sorted_q(1:6))+1,-sorted_ratio_diff(1:6),sorted_TF_array(1:6));

xlabel('-log10(p)')
ylabel('mean ratio - mean ratio')


disp(length(find(q < 0.01)))


%% Save supplementary data 1
P_value = q';
Ratio_difference = -ratio_diff';
T = table(TF_array,P_value,Ratio_difference);
writetable(T,'supplementary_data_1_2.xls','FileType','spreadsheet')


%% **************************** INDIVIDUAL TFS ***************************
%% Find positions for specific transcription factor based on motif sequence
TF_name = 'Reb1'; % Abf1, Rap1, Reb1
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


%% Produce sequence logos for all sites
motif_sequences = cell(1,length(site.chr));
for site_index = 1:length(site.chr)
    if strcmp(site.strand(site_index),'+')
        motif_sequences{site_index} = char(site.motif(site_index));
    else
        motif_sequences{site_index} = seqrcomplement(char(site.motif(site_index)));
    end
end
seqlogo(motif_sequences);


%% Calculate ratio around each predicted TF site
site.ratio = compute_property_around_sites(yeast_genome,cell_naked_ratio,22,site);


%% Investigate if motif sites overlap with DNase I footprints
site.footprint = check_genomic_overlap({yeast_genome.Header},site,dnase_footprints);


%% Find overlap with ChIP-peaks
if strcmp(TF_name,'Reb1')
    chip.chr = Reb1_chip.Var1;
    chip.peak = Reb1_chip.Var2;
    chip.width = Reb1_chip.Var4;
    chip.start = chip.peak - ceil(chip.width / 2);
    chip.stop = chip.peak + ceil(chip.width / 2);

    site.chip = check_genomic_overlap({yeast_genome.Header},site,chip);
end


%% Find distance to, and gene expression of, closest downstream genes
site.distance_to_gene = nan(1,length(site.chr));
for site_index = 1:length(site.chr)
    site_center = (site.start(site_index) + site.stop(site_index)) / 2;
    gene_strand_current_chromosome = genes.strand(strcmp(string(genes.chr),site.chr(site_index)));
    gene_start_current_chromosome = genes.start(strcmp(string(genes.chr),site.chr(site_index)));
    
    index = find((gene_strand_current_chromosome == '+' & gene_start_current_chromosome > site_center) | (gene_strand_current_chromosome == '-' & gene_start_current_chromosome < site_center));
    distances = abs(gene_start_current_chromosome(index) - site_center);
    site.distance_to_gene(site_index) = min(distances);
end


%% Re-format the DNase I cleavage data
cleavage_new = cell(1,length(yeast_genome));
for chromosome_index = 1:length(yeast_genome)
    cleavage_new{chromosome_index} = zeros(1,length(yeast_genome(chromosome_index).Sequence));

    positions_current_chromosome = dnase_tagcounts.start(strcmp(dnase_tagcounts.chr,{yeast_genome(chromosome_index).Header}));
    cleavage_current_chromosome = dnase_tagcounts.cleavage(strcmp(dnase_tagcounts.chr,{yeast_genome(chromosome_index).Header}));

    cleavage_new{chromosome_index}(positions_current_chromosome+1) = cleavage_current_chromosome;
end


%% Calculate DNase I cleavage around each predicted TF site
site.dnase = compute_property_around_sites(yeast_genome,cleavage_new,22,site);


%%
disp(min(min(site.ratio(23:23+length(char(site.motif(1)))-1,find(site.footprint)))))


%% Sort sites based on ratio
site.site_ratio = mean(abs(site.ratio(23:23+length(char(site.motif(1)))-1,:)),1,'omitnan');
[~,sorting_index] = sort(site.site_ratio,'descend','MissingPlacement','last');


fields = fieldnames(site);
for i = 1:length(fieldnames(site))
    site.(fields{i}) = site.(fields{i})(:,sorting_index);
end


%% Plot mean signal of the motif sites
idx = site.footprint;
motif_length = length(char(site.motif(1)));

figure()
ratio_cluster1 = mean(site.ratio(:,idx == 1),2,'omitnan');
ratio_cluster2 = mean(site.ratio(:,idx == 0),2,'omitnan');
ratio_mean = [ratio_cluster1'; ratio_cluster2'];
xpos = 1:size(site.ratio,1);
bar(xpos,ratio_mean,BarWidth=1)
side_width = (25-motif_length)/2;
xlim([23-floor(side_width)-0.5 23+motif_length-1+ceil(side_width)+0.5]);
xticks(xpos)
xlabel('Position relative binding motif (b)',Fontsize=10);
ylabel('Mean log_2(f_{sync}/f_{naked})',Fontsize=10);
legend(['DNase I overlap, n=' num2str(length(find(idx == 1)))],['No overlap, n=' num2str(length(find(idx == 0)))],fontname='helvetica',fontsize=14)
ylim([-4 4]);


%% Plot heatmap
color_map = readmatrix('../colormaps/diverging_colormap_redblue_new2','FileType','text');

figure()
subplot(1,7,[1 4])
imAlpha=ones(size(site.ratio'));
imAlpha(isnan(site.ratio')) = 0;
image(site.ratio',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.8 0.8 0.8]);
colormap(flip(color_map));
clim([-3 3]);
colorbar(Location='westoutside');
sgtitle('Log_2(f_{sync}/f_{naked}) for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:25:51;
ax.XTickLabel = -25:25:25;
xlim([23-22-0.5 23+motif_length-1+22+0.5])

subplot(1,7,5)
barh(flip(site.footprint),Facecolor='k');
xlabel('Footprints')

if strcmp(TF_name,'Reb1')
    subplot(1,7,6)
    barh(flip(site.chip),Facecolor='k');
    xlabel('ChIP')
end

subplot(1,7,7)
barh(flip(site.distance_to_gene),Facecolor='k')
xlim([0 5000])


%% Plot DNase I cleavage heatmap    
figure()
image(site.dnase',CDataMapping='scaled');
c = gray;
colormap(flip(c));
clim([0 8]);
colorbar(Location='westoutside');
sgtitle('DNase I cleavage for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:5:51;
ax.XTickLabel = -25:5:25;


%% Compare relationship between Deam-seq ratio and DNase footprints / ChIP-exo peaks
p_dnase = ranksum(site.site_ratio(logical(site.footprint)),site.site_ratio(~logical(site.footprint)));
p_dnase = p_dnase * length(site.site_ratio);
fprintf('p dnase: %d \n',p_dnase)

if strcmp(TF_name,'Reb1')
    p_chip = ranksum(site.site_ratio(logical(site.chip)),site.site_ratio(~logical(site.chip)));
    p_chip = p_chip * length(site.site_ratio);
    fprintf('p chip: %d \n',p_chip)
end

% For DNase:
median_ratio_footprint = median(site.site_ratio(logical(site.footprint)),'omitnan');
median_ratio_no_footprint = median(site.site_ratio(~logical(site.footprint)),'omitnan');
iqr_ratio_footprint = iqr(site.site_ratio(logical(site.footprint)));
iqr_ratio_no_footprint = iqr(site.site_ratio(~logical(site.footprint)));

% For ChIP:
if strcmp(TF_name,'Reb1')
    median_ratio_chip = median(site.site_ratio(logical(site.chip)),'omitnan');
    median_ratio_no_chip = median(site.site_ratio(~logical(site.chip)),'omitnan');
    iqr_ratio_chip = iqr(site.site_ratio(logical(site.chip)));
    iqr_ratio_no_chip = iqr(site.site_ratio(~logical(site.chip)));
end


% CDF-plots:
figure()
cdfplot(site.site_ratio(logical(site.footprint)))
hold on
cdfplot(site.site_ratio(~logical(site.footprint)))
legend('Overlap','No overlap','Location','northwest')

if strcmp(TF_name,'Reb1')
    figure()
    cdfplot(site.site_ratio(logical(site.chip)))
    hold on
    cdfplot(site.site_ratio(~logical(site.chip)))
    legend('Overlap','No overlap','Location','northwest')
end
