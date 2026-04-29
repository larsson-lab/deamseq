addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat

cell_naked_ratio = mutation_ratio{1}.penta_norm;

% Load gene expression data:
expression = readtable('../data/expression/FPKM_NR_CR.txt','FileType','text','ReadVariableNames',true,'Delimiter','\t');


%% Calculate ratio relative TSSs
genes.stop = genes.start;

% Compute ratio around TSS-1000 and TSS+1000:
ratio_gene = compute_property_around_sites(yeast_genome,cell_naked_ratio,1000,genes)';


%% Sort the ratios relative TSSs according to gene expression, and remove rows for genes with no expression
% Map expression data to genes:
[~,gene_data_index,expression_data_index] = intersect(genes.name,expression.gene_short_name);

% Expression level of each gene:
gene_expression_level = nan(1,length(genes.start));
gene_expression_level(gene_data_index) = (expression.s_1_NR_A_FPKM(expression_data_index) + expression.s_2_NR_B_FPKM(expression_data_index)) / 2;

% Sort genes according to expression in descending order:
[gene_expression_level_sorted,index] = sort(gene_expression_level,'descend','MissingPlacement','last');

% Sort ratios according to gene expression, and remove rows for genes with no expression:
ratio_sorted = ratio_gene(index,:);
index_no_expression = isnan(gene_expression_level_sorted);
ratio_sorted(index_no_expression,:) = [];
gene_expression_level_sorted(:,index_no_expression) = [];


%% Sort the positions and genes into bins of 20 each:
number_xbins = ceil(2001/20);
number_ybins = ceil(length(ratio_sorted)/20);
bin_values_max_abs = zeros(number_xbins,number_ybins);
current_value_max_abs = -inf(number_xbins,number_ybins);
expression_sum = zeros(1,number_ybins);
for gene_index = 1:length(ratio_sorted)
    bin_ynr = ceil(gene_index/20);
    for position = 1:2001
        bin_xnr = ceil(position/20);

        bin_values_max_abs(bin_xnr,bin_ynr) = max(abs(ratio_sorted(gene_index,position)),current_value_max_abs(bin_xnr,bin_ynr));
        current_value_max_abs(bin_xnr,bin_ynr) = bin_values_max_abs(bin_xnr,bin_ynr);
    end
    expression_sum(bin_ynr) = expression_sum(bin_ynr) + gene_expression_level_sorted(gene_index);
end
average_expression_per_bin = expression_sum ./ number_ybins;


%% Produce heatmaps
color_map = readmatrix('../colormaps/greenblue_colormap_new3.txt','FileType','text');

bin_values_max_abs(bin_values_max_abs==-Inf) = nan;

figure()
subplot(1,5,[1 4])
imAlpha=ones(size(bin_values_max_abs'));
imAlpha(isnan(bin_values_max_abs')) = 0;
image(bin_values_max_abs',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[0.59 0.14 0.75]);
colormap(color_map);
clim([0 4]);
colorbarTickLabels = {'0','0.5','1','1.5','2','2.5','3','3.5',[char(8805) '4']};
colorbar(Location='westoutside',TickLabels=colorbarTickLabels);
sgtitle('Maximum absolute ratio of log_2(f_{sync}/f_{naked})',Fontsize=14);
ax = gca();
ax.XTick = 1:10:101;
ax.XTickLabel = -1000:200:1000;

subplot(1,5,5)
barh(flip(average_expression_per_bin),Facecolor=[0.5 0.5 0.5]);
set(gca,'xscale','log')
xlim([0 100000]);
set(gca,'YTickLabel',[]);


%% Bin the whole genome into 20 bp windows
% Extract maximum abs. signal of each bin, and exact genomic position of that signal:
binned_signal = cell(1,length(yeast_genome));
binned_signal_pos = cell(1,length(yeast_genome));
for chromosome_index = 1:length(yeast_genome)
    edited_ratio = [cell_naked_ratio{chromosome_index} nan(1,20-mod(length(cell_naked_ratio{chromosome_index}),20))];
    nr_bins = length(edited_ratio)/20;
    ratio_matrix = reshape(edited_ratio,[20,nr_bins]);

    for bin_index = 1:nr_bins
        binned_signal{chromosome_index}(bin_index) = max(abs(ratio_matrix(:,bin_index)));
        index_max = abs(ratio_matrix(:,bin_index)) == max(abs(ratio_matrix(:,bin_index)));

        index = nan;
        if isscalar(find(index_max))
            index = 20 * (bin_index - 1) + find(index_max) - 1;
        end

        binned_signal_pos{chromosome_index}(bin_index) = index;
    end
end


%% Iterate through the binned data and extract the strongest signals
signal_cutoffs = [1 2 3 4 inf];
dnase_counts = zeros(length(signal_cutoffs),401);
for cutoff_index = 1:length(signal_cutoffs)-1
    for chromosome_index = 1:length(yeast_genome)
        strongest_pos = binned_signal_pos{chromosome_index}(binned_signal{chromosome_index} >= signal_cutoffs(cutoff_index) & binned_signal{chromosome_index} < signal_cutoffs(cutoff_index+1));
    
        footprint_start_current_chr = dnase_footprints.start(strcmp({yeast_genome(chromosome_index).Header},dnase_footprints.chr));
        footprint_stop_current_chr = dnase_footprints.stop(strcmp({yeast_genome(chromosome_index).Header},dnase_footprints.chr));
        footprint_boolean = zeros(1,length(yeast_genome(chromosome_index).Sequence));
        for footprint_index = 1:length(footprint_start_current_chr)
            footprint_boolean(footprint_start_current_chr(footprint_index):footprint_stop_current_chr(footprint_index)) = 1;
        end
    
        for i = 1:length(strongest_pos)
            if strongest_pos(i) <= 200 || strongest_pos(i) + 200 >= length(footprint_boolean)
                continue
            end
            dnase_counts(cutoff_index,:) = dnase_counts(cutoff_index,:) + footprint_boolean(strongest_pos(i)-200:strongest_pos(i)+200);
        end
    end
end


%% Plot the normalized counts for each signal strength
figure()
hold on
for cutoff_index = 1:length(signal_cutoffs)-1
    plot(dnase_counts(cutoff_index,:) / sum(dnase_counts(cutoff_index,:)))
end
xticks(1:100:401);
xticklabels(-200:100:200);
xlim([1 401]);
ylim([0 20e-3])
xlabel('Position relative strong deam-seq peak (bp)');
ylabel('Normalized DNase footprint counts');
legend('1-2','2-3','3-4','>=4')


%% Find genome-wide strongest positive and negative signals
current_max = -Inf;
current_min = Inf;
for chromosome_index = 1:length(yeast_genome)
    max_signal = max(current_max,max(cell_naked_ratio{chromosome_index}));
    current_max = max_signal;

    min_signal = min(current_min,min(cell_naked_ratio{chromosome_index}));
    current_min = min_signal;
end


% Specific region:
min(mutation_ratio{1}.penta_norm{strcmp({yeast_genome.Header},'chrVIII')}(350000:400000))