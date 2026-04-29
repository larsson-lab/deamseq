addpath('../functions');

%% Load data
load '../workspace.mat'
load ../extra_data.mat

% Interested in sync/naked ratio:
cell_naked_ratio = mutation_ratio{1}.penta_norm;

% Load nucleosome data:
nucleosomes = readtable('../data/nucleosomes/lift_over/41586_2012_BFnature11142_MOESM263_ESM.reformatted.bed','FileType','text','ReadVariableNames',false,'Delimiter','\t');
nucleosomes = removevars(nucleosomes,"Var3");
nucleosomes = renamevars(nucleosomes,["Var1","Var2","Var5"],["chr","pos","score"]);


%% Calculate ratio relative TSSs
genes.stop = genes.start;

% Compute ratio around TSS-1000 and TSS+1000:
ratio_gene = compute_property_around_sites(yeast_genome,cell_naked_ratio,1000,genes)';


%% Plot mean ratio relative to downstream gene TSSs
figure()
subplot(2,1,1)
xpos = -1000:1000;
ratio_gene_mean = mean(ratio_gene,'omitnan');
plot(xpos,ratio_gene_mean,Color=[0.28 0.54 0.37],LineWidth=1);
xlabel('Position relative TSS',Fontsize=10);
ylabel('Mean log2(f_{sync}/f_{naked})',Fontsize=10);
ylim([-0.2 0.1]);


subplot(2,1,2)
xpos = -200:200;
plot(xpos,ratio_gene_mean(801:1201),Color=[0.26 0.54 0.27],LineWidth=1);
xlabel('Position relative TSS',Fontsize=10);
ylabel('Mean log2(f_{sync}/f_{naked})',Fontsize=10);
ylim([-0.2 0.1]);


%% Plot the FFT of the mean ratio relative downstream gene TSSs
ratio_array = {ratio_gene,ratio_gene(:,1:1000),ratio_gene(:,1002:2001)};
title_array = {'Whole region','Upstream of TSS','Downstream of TSS'};


for plot_index = 1:length(ratio_array)
    fourier_ratio = fft(mean(ratio_array{plot_index},'omitnan'));
    number_data_points = size(fourier_ratio,2);
    
    figure()
    xpos = 0:floor(number_data_points/2)-1;
    signal = abs(fourier_ratio(1:floor(number_data_points/2)));
    plot(xpos,signal,Color='k',LineWidth=1);
    axis([0 floor(number_data_points/2) 0 22]);
    hold on
    
    period = number_data_points ./ xpos;
    
    % Peak at period ~ 150:
    peak1 = max(signal(period >= 100 & period <= 170));
    index1 = find(signal == peak1);
    if peak1 >= 1
        line([index1 index1+30],[peak1+0.2 peak1+0.2],Color='k');
        text(index1+33,peak1+0.2,string(period(index1)),Color='k',Fontsize=14);
    end
    % Peak at period ~ 10:
    peak2 = max(signal(period <= 15 & period >= 5));
    index2 = find(signal == peak2);
    if peak2 >= 1
        line([index2 index2+30],[peak2+0.2 peak2+0.2],Color='k');
        text(index2+33,peak2+0.2,string(period(index2)),Color='k',Fontsize=14);
    end
    
    ticks = zeros(1,floor(number_data_points/2)/100);
    for i = 1:floor(number_data_points/2)/100
        ticks(i) = period(i*floor(number_data_points/2)/10 - (floor(number_data_points/2)/10 - 1));
    end
    set(gca,'XTick',0:floor(number_data_points/2)/10:length(period)-1,'XTickLabel',ticks(:));
    xlabel('Period (b)',Fontsize=10);
    ylabel('FFT[mean log_2(f_{sync}/f_{naked})]',Fontsize=10);
    sgtitle(title_array{plot_index},Fontsize=14);
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


ratio_nucleosome = compute_property_around_sites(yeast_genome,cell_naked_ratio,1000,nucleosomes);


%% Plot mean ratio relative nucleosome positions for top and bottom quintiles of nucleosome scores
ratio_nucleosome_highocc = ratio_nucleosome(:,1:round(length(ratio_nucleosome)/5));
ratio_nucleosome_lowocc = ratio_nucleosome(:,round(4*length(ratio_nucleosome)/5):end);

ratio_nucleosome_highocc_mean = mean(ratio_nucleosome_highocc,2,'omitnan');
ratio_nucleosome_lowocc_mean = mean(ratio_nucleosome_lowocc,2,'omitnan');
figure()
subplot(2,1,1)

xpos = -1000:1000;
plot(xpos,ratio_nucleosome_lowocc_mean,Color=[0 0.29 0.44],LineWidth=1);
hold on
plot(xpos,ratio_nucleosome_highocc_mean,Color=[0.51 0.65 0.28],LineWidth=1);
xlabel('Position relative nucleosome',Fontsize=10);
ylabel('Mean log_2(f_{sync}/f_{naked})',Fontsize=10);
legend('5th quintile nucleosome occupancy','1st quintile nucleosome occupancy',fontsize=7,fontname='helvetica',edgecolor='none');
axis([-1000 1000 -0.25 0.25]);

subplot(2,1,2)
plot(xpos,ratio_nucleosome_lowocc_mean,Color=[0 0.29 0.44],LineWidth=1);
hold on
plot(xpos,ratio_nucleosome_highocc_mean,Color=[0.51 0.65 0.28],LineWidth=1);
outwardPositions = -210:10.5:210;
xline(outwardPositions,Color=[0.5 0.5 0.5],LineStyle=":");
xlabel('Position relative nucleosome',Fontsize=10);
ylabel('Mean log_2(f_{sync}/f_{naked})',Fontsize=10);
legend('5th quintile nucleosome occupancy','1st quintile nucleosome occupancy',fontsize=7,fontname='helvetica',edgecolor='none');
axis([-200 200 -0.25 0.25]);


%% Plot the FFT of the mean ratio relative (1st quintile occupancy) nucleosome positions
fourier_ratio = fft(ratio_nucleosome_highocc_mean);
number_data_points = size(fourier_ratio,1);

figure()
xpos = 0:floor(number_data_points/2)-1;
signal = abs(fourier_ratio(1:floor(number_data_points/2)));
plot(xpos,signal,Color='k',LineWidth=1);
axis([0 floor(number_data_points/2) 0 35]);
xlim([0 floor(number_data_points/2)]);
hold on

period = number_data_points ./ xpos;

% Peak at period ~ 150:
peak1 = max(signal(period >= 100 & period <= 170));
index1 = find(signal == peak1);
if peak1 >= 1
    line([index1 index1+30],[peak1+0.05 peak1+0.05],Color='k');
    text(index1+33,peak1+0.05,string(period(index1)),Color='k',Fontsize=14);
end

% Peak at period ~ 10:
peak2 = max(signal(period <= 15 & period >= 5));
index2 = find(signal == peak2);
if peak2 >= 1
    line([index2 index2+30],[peak2+0.05 peak2+0.05],Color='k');
    text(index2+33,peak2+0.05,string(period(index2)),Color='k',Fontsize=14);
end

ticks = zeros(1,floor(number_data_points/2)/100);
for i = 1:floor(number_data_points/2)/100
    ticks(i) = period(i*floor(number_data_points/2)/10 - (floor(number_data_points/2)/10 - 1));
end
set(gca,'XTick',0:floor(number_data_points/2)/10:length(period)-1,'XTickLabel',ticks(:));
xlabel('Period (bp)',Fontsize=10);
ylabel('FFT[mean log_2(f_{sync}/f_{naked})]',Fontsize=10);


%% Sort the nucleosomes into bins of 200 nucleosomes and 20 positions each:
number_xbins = ceil(2001/20);
number_ybins = ceil(length(nucleosomes.chr)/200);
bin_values_sum = nan(number_xbins,number_ybins);
number_bins = zeros(number_xbins,number_ybins);
occupancy_sum = zeros(1,number_ybins);
for nucleosome_index = 1:length(nucleosomes.chr)
    bin_ynr = ceil(nucleosome_index/200);
    for position = 1:2001
        bin_xnr = ceil(position/20);

        % Sum all non-NaN ratios in the bin:
        if isnan(ratio_nucleosome(position,nucleosome_index))
            continue
        end
        if isnan(bin_values_sum(bin_xnr,bin_ynr))
            bin_values_sum(bin_xnr,bin_ynr) = ratio_nucleosome(position,nucleosome_index);
        else
            bin_values_sum(bin_xnr,bin_ynr) = bin_values_sum(bin_xnr,bin_ynr) + ratio_nucleosome(position,nucleosome_index);
        end
        number_bins(bin_xnr,bin_ynr) = number_bins(bin_xnr,bin_ynr) + 1;
    end

    % Sum all nucleosome occupancies for each y-bin:
    occupancy_sum(bin_ynr) = occupancy_sum(bin_ynr) + nucleosomes.score(nucleosome_index);
end
bin_values_mean = bin_values_sum ./ number_bins;
average_occupancy_per_bin = occupancy_sum ./ number_ybins;


%% Produce heatmap with mean ratios
color_map = readmatrix('../colormaps/diverging_colormap_redblue_new2.txt','FileType','text');

figure()
subplot(1,5,[1 4])
imAlpha=ones(size(bin_values_mean'));
imAlpha(isnan(bin_values_mean')) = 0;
image(bin_values_mean',CDataMapping='scaled',AlphaData=imAlpha);
set(gca,'Color',[1 1 1]);

colormap(flip(color_map));
clim([-0.1 0.1]);
colorbar(Location='westoutside');
sgtitle('Mean log_2(f_{sync}/f_{naked}) per 200 nucl. and 20 pos.',Fontsize=14);
ax = gca();
ax.XTick = 1:10:201;
ax.XTickLabel = -1000:200:1000;

subplot(1,5,5)
barh(flip(average_occupancy_per_bin),Facecolor=[0.5 0.5 0.5]);
set(gca,'xscale','log')
xlim([0 10]);
set(gca,'YTickLabel',[]);
