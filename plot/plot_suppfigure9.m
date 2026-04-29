addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat

% Load gene expression data:
expression = readtable('../data/expression/FPKM_NR_CR.txt','FileType','text','ReadVariableNames',true,'Delimiter','\t');

% Load binding data for sites:
Reb1_chip_depth = readtable('../data/ChIP/Reb1_Rap1/SRR346400_depth.tsv','FileType','text','ReadVariableNames',false,'Delimiter','\t');

% Load motif sequences for transcription factors:
motifs = readtable('../data/motifs/TF_motif_sequences_short.txt','FileType','text','ReadVariableNames',false,'Delimiter',',');
motifs = renamevars(motifs,["Var1","Var2","Var3","Var4"],["Name","Seq","Sym_seq","Sym_seq_rcomp"]);

% Load nucleosome data:
nucleosomes = readtable('../data/nucleosomes/lift_over/sacCer3.41586_2012_BFnature11142_MOESM263_ESM.reformatted.bed','FileType','text','ReadVariableNames',false,'Delimiter','\t');


%% Find positions for Reb1 based on motif sequence
TF_name = 'Reb1';
TF_motif_index = find(strcmp(motifs.Name,TF_name));
motif_sequence = motifs.Seq{TF_motif_index};

% Find motif positions:
site = find_seq_motif_sites(yeast_genome,motif_sequence,motifs.Sym_seq{TF_motif_index});

fprintf(['Number of motif sites for ' TF_name ': %d\n'],length(site.chr));


%% Remove all sites that are not included in the whitelisted regions
outside_whitelisted = false(1,length(site.chr));
for site_idx = 1:length(site.chr)
    current_chr = strcmp(site.chr(site_idx),{yeast_genome.Header});
    outside_whitelisted(site_idx) = any(idx_cov_ok{current_chr}(site.start(site_idx):site.stop(site_idx)) == 0);
end
site.chr = site.chr(~outside_whitelisted);
site.start = site.start(~outside_whitelisted);
site.stop = site.stop(~outside_whitelisted);
site.strand = site.strand(~outside_whitelisted);
site.motif = site.motif(~outside_whitelisted);


%% Investigate if motif sites overlap with DNase I footprints
site.footprint = check_genomic_overlap({yeast_genome.Header},site,dnase_footprints);


%% Calculate gene expression of closest downstream TSS for each gene
% Map expression data to genes:
[~,gene_data_index,expression_data_index] = intersect(genes.name,expression.gene_short_name);

% Expression level of each gene:
gene_expression_level = nan(1,length(genes.start));
gene_expression_level(gene_data_index) = (expression.s_1_NR_A_FPKM(expression_data_index) + expression.s_2_NR_B_FPKM(expression_data_index)) / 2;


site.expr = nan(1,length(site.chr));
for site_idx = 1:length(site.chr)
    genes_idx = genes.chr == site.chr(site_idx);
    genes_strand = genes.strand(genes_idx);
    genes_start = genes.start(genes_idx);

    % Find closest downstream gene:
    if strcmp(site.strand(site_idx),'+')
        closest_gene = min(genes_start(genes_start >= site.start(site_idx)));
    elseif strcmp(site.strand(site_idx),'-')
        closest_gene = max(genes_start(genes_start <= site.start(site_idx)));
    end

    if isempty(closest_gene)
        continue
    end

    closest_gene_idx = find(genes_start == double(closest_gene));
    site.expr(site_idx) = gene_expression_level(closest_gene_idx(1));
end


%% Re-format the DNase I cleavage data
cleavage_new = cell(1,length(yeast_genome));
for chr_idx = 1:length(yeast_genome)
    cleavage_new{chr_idx} = zeros(1,length(yeast_genome(chr_idx).Sequence));

    positions_current_chromosome = dnase_tagcounts.start(strcmp(dnase_tagcounts.chr,{yeast_genome(chr_idx).Header}));
    cleavage_current_chromosome = dnase_tagcounts.cleavage(strcmp(dnase_tagcounts.chr,{yeast_genome(chr_idx).Header}));

    cleavage_new{chr_idx}(positions_current_chromosome+1) = cleavage_current_chromosome;
end


%% Calculate DNase I cleavage around each predicted TF site
site.dnase = compute_property_around_sites(yeast_genome,cleavage_new,1000,site);


%% Calculate ratios (Reb1aa and Reb1aa+rap) around each predicted TF site
site.ratio_aa = compute_property_around_sites(yeast_genome,mutation_ratio{4}.penta_norm,22,site);
site.ratio_aa_rap = compute_property_around_sites(yeast_genome,mutation_ratio{5}.penta_norm,22,site);


%% Re-format the nucleosome data
nucleosome_chrom_pos = cell(1,length(yeast_genome));
for chr_idx = 1:length(yeast_genome)
    nucleosome_chrom_pos{chr_idx} = zeros(1,length(yeast_genome(chr_idx).Sequence));

    nucleosome_current_chr = nucleosomes.Var2(strcmp(nucleosomes.Var1,string(yeast_genome(chr_idx).Header)));
    nucleosome_score_current_chr = nucleosomes.Var5(strcmp(nucleosomes.Var1,string(yeast_genome(chr_idx).Header)));
    nucleosome_chrom_pos{chr_idx}(nucleosome_current_chr) = nucleosome_score_current_chr;
end


%% Calculate nucleosome occupancy around each predicted TF site
site.nucleosome = compute_property_around_sites(yeast_genome,nucleosome_chrom_pos,1000,site);


%% Re-format the ChIP-exo coverage data
chip_coverage = cell(1,length(yeast_genome));
for chr_idx = 1:length(yeast_genome)
    chip_coverage{chr_idx} = zeros(1,length(yeast_genome(chr_idx).Sequence));
    coverage_current_chromosome = Reb1_chip_depth.Var3(strcmp(string(Reb1_chip_depth.Var1),{yeast_genome(chr_idx).Header}));
    chip_coverage{chr_idx}(1:length(yeast_genome(chr_idx).Sequence)) = coverage_current_chromosome;
end


%% Calculate ChIP-exo coverage around each predicted TF site
site.chip_cov = compute_property_around_sites(yeast_genome,chip_coverage,1000,site);


%% Sort sites based on ratio
ratio_motif_aa = mean(abs(site.ratio_aa(23:23+length(char(site.motif(1)))-1,:)),'omitnan');
ratio_motif_aa_rap = mean(abs(site.ratio_aa_rap(23:23+length(char(site.motif(1)))-1,:)),'omitnan');
site.ratio_motif_diff = ratio_motif_aa - ratio_motif_aa_rap;

[~,sorting_index] = sort(site.ratio_motif_diff,'descend','MissingPlacement','last');


fields = fieldnames(site);
for i = 1:length(fieldnames(site))
    site.(fields{i}) = site.(fields{i})(:,sorting_index);
end


%% Find DNase I unsupported sites with strong Deam-seq signal diff
idx_not_supported = find(~site.footprint);
idx_not_supported = idx_not_supported(mean(site.dnase(701:1301,idx_not_supported),1,'omitnan') <= 5);
idx_strong_signal = find(site.ratio_motif_diff >= 0.75);
idx = intersect(idx_strong_signal,idx_not_supported);
idx2 = idx_strong_signal(~ismember(idx_strong_signal,idx));


%% Plot heatmap for Reb1aa (DNase I supported and unsupported sites)
motif_length = length(char(site.motif(1)));
color_map = readmatrix('../colormaps/diverging_colormap_redblue_new2','FileType','text');

% Non-supported:
figure()
subplot(1,4,[1 2])
imAlpha=ones(size(site.ratio_aa(:,idx)'));
imAlpha(isnan(site.ratio_aa(:,idx)')) = 0;
image(site.ratio_aa(:,idx)',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.8 0.8 0.8]);
colormap(flip(color_map));
clim([-3 3]);
colorbar(Location='westoutside');
sgtitle('Log_2(f_{sync}/f_{naked}) for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:25:51;
ax.XTickLabel = -25:25:25;
xlim([23-22-0.5 23+motif_length-1+22+0.5])

subplot(1,4,[3 4])
imAlpha=ones(size(site.ratio_aa_rap(:,idx)'));
imAlpha(isnan(site.ratio_aa_rap(:,idx)')) = 0;
image(site.ratio_aa_rap(:,idx)',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.8 0.8 0.8]);
colormap(flip(color_map));
clim([-3 3]);
colorbar(Location='westoutside');
sgtitle('Log_2(f_{sync}/f_{naked}) for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:25:51;
ax.XTickLabel = -25:25:25;
xlim([23-22-0.5 23+motif_length-1+22+0.5])


% Supported:
figure()
subplot(1,4,[1 2])
imAlpha=ones(size(site.ratio_aa(:,idx2)'));
imAlpha(isnan(site.ratio_aa(:,idx2)')) = 0;
image(site.ratio_aa(:,idx2)',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.8 0.8 0.8]);
colormap(flip(color_map));
clim([-3 3]);
colorbar(Location='westoutside');
sgtitle('Log_2(f_{sync}/f_{naked}) for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:25:51;
ax.XTickLabel = -25:25:25;
xlim([23-22-0.5 23+motif_length-1+22+0.5])

subplot(1,4,[3 4])
imAlpha=ones(size(site.ratio_aa_rap(:,idx2)'));
imAlpha(isnan(site.ratio_aa_rap(:,idx2)')) = 0;
image(site.ratio_aa_rap(:,idx2)',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.8 0.8 0.8]);
colormap(flip(color_map));
clim([-3 3]);
colorbar(Location='westoutside');
sgtitle('Log_2(f_{sync}/f_{naked}) for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:25:51;
ax.XTickLabel = -25:25:25;
xlim([23-22-0.5 23+motif_length-1+22+0.5])


%% Plot sequence logos for the two sets
sequence_nodnase = cell(1,length(idx));
for i = 1:length(idx)
    current_chr = strcmp(site.chr(idx(i)),{yeast_genome.Header});
    if strcmp(site.strand(idx(i)),'+')
        sequence_nodnase{i} = yeast_genome(current_chr).Sequence(site.start(idx(i))-5:site.stop(idx(i))+5);
    elseif strcmp(site.strand(idx(i)),'-')
        sequence_nodnase{i} = seqrcomplement(yeast_genome(current_chr).Sequence(site.start(idx(i))-5:site.stop(idx(i))+5));
    end
end

sequence_dnase = cell(1,length(idx2));
for i = 1:length(idx2)
    current_chr = strcmp(site.chr(idx2(i)),{yeast_genome.Header});
    if strcmp(site.strand(idx2(i)),'+')
        sequence_dnase{i} = yeast_genome(current_chr).Sequence(site.start(idx2(i))-5:site.stop(idx2(i))+5);
    elseif strcmp(site.strand(idx2(i)),'-')
        sequence_dnase{i} = seqrcomplement(yeast_genome(current_chr).Sequence(site.start(idx2(i))-5:site.stop(idx2(i))+5));
    end
end

seqlogo(sequence_nodnase)
seqlogo(sequence_dnase)


%% Plot the DNase cleavage data
% Non-supported:
figure()
subplot(1,5,[1 2])
image(site.dnase(:,idx)',CDataMapping='scaled');
c = gray;
colormap(flip(c));
clim([0 8]);
colorbar(Location='westoutside');
sgtitle('DNase I cleavage for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
xlim([700 1301])

subplot(1,5,[3 4])
plot(mean(site.dnase(:,idx)'));
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
axis([700 1301 0 10])

% Supported:
figure()
subplot(1,5,[1 2])
image(site.dnase(:,idx2)',CDataMapping='scaled');
c = gray;
colormap(flip(c));
clim([0 8]);
colorbar(Location='westoutside');
sgtitle('DNase I cleavage for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
xlim([700 1301])

subplot(1,5,[3 4])
plot(mean(site.dnase(:,idx2)','omitnan'));
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
axis([700 1301 0 10])


%% Look at nucleosomes relative the predicted Reb1 sites
% Non-supported:
figure()
plot(mean(site.nucleosome(:,idx)'));
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
axis([700 1301 0 0.4])
title('Nucleosome occupancy for predicted sites',Fontsize=14);

% Supported
figure()
plot(mean(site.nucleosome(:,idx2)','omitnan'));
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
axis([700 1301 0 0.4])
title('Nucleosome occupancy for predicted sites',Fontsize=14);


%% Plot the coverage for the ChIP-exo data
figure()
image(site.chip_cov(:,idx)',CDataMapping='scaled');
c = gray;
colormap(flip(c));
clim([0 100]);
colorbar(Location='westoutside');
sgtitle('ChIP-exo coverage for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
xlim([700 1301])

figure()
image(site.chip_cov(:,idx2)',CDataMapping='scaled');
c = gray;
colormap(flip(c));
clim([0 100]);
colorbar(Location='westoutside');
sgtitle('ChIP-exo coverage for predicted sites',Fontsize=14);
ax = gca();
ax.XTick = 1:100:2001;
ax.XTickLabel = -1000:100:1000;
xlim([700 1301])


%% Plot distribution of maximum DNase I cleavages
p = ranksum(mean(site.dnase(700:1301,idx)),mean(site.dnase(700:1301,idx2)));
disp(p)

figure()
cdfplot(mean(site.dnase(700:1301,idx)))
hold on
cdfplot(mean(site.dnase(700:1301,idx2)))
xlabel('Mean DNase I cleavage')
legend('DNase I non-supported','DNase I supported','Location','southeast')


%% Plot distribution of mean nucleosome occupancies per window
p = ranksum(mean(site.nucleosome(700:1301,idx),'omitnan'),mean(site.nucleosome(700:1301,idx2),'omitnan'));
disp(p)

figure()
cdfplot(mean(site.nucleosome(700:1301,idx),'omitnan'))
hold on
cdfplot(mean(site.nucleosome(700:1301,idx2),'omitnan'))
xlabel('Sum nucleosome occupancy')
legend('DNase I non-supported','DNase I supported','Location','southeast')

