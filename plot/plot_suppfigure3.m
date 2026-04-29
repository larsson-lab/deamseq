addpath('../functions');

%% Load data
load '../workspace.mat'
load ../extra_data.mat


%% Plot distribution of ratios
figure()
subplot(2,2,[1 2])
hold on
for plot_index = [1 2 3 6]
    whole_ratio = [];
    for chromosome_index = 1:length(yeast_genome)
        whole_ratio = [whole_ratio mutation_ratio{plot_index}.penta_norm{chromosome_index}];
    end
    
    ratio_sorted = sort(whole_ratio,'ascend','MissingPlacement','last');
    plot(ratio_sorted,'LineWidth',2);
    if plot_index == 1
        disp(min(ratio_sorted))
    end
end
yline(0,'LineWidth',1)
xlim([0 length(find(~isnan(ratio_sorted)))])
ylim([-7 7])
ylabel('Cell/naked deam-seq ratio (log2)')
legend('Cell/naked pooled','Cell1/naked1','Cell2/naked2','Cell1/cell2')


subplot(2,2,3)
hold on
for plot_index = [1 2 3 6]
    whole_ratio = [];
    for chromosome_index = 1:length(yeast_genome)
        whole_ratio = [whole_ratio mutation_ratio{plot_index}.penta_norm{chromosome_index}];
    end
    
    ratio_sorted = sort(whole_ratio,'ascend','MissingPlacement','last');
    plot(ratio_sorted,'LineWidth',2);
end
yline(0,'LineWidth',1)
xlim([0 10000])
ylim([-7 7])
ylabel('Cell/naked deam-seq ratio (log2)')
legend('Cell/naked pooled','Cell1/naked1','Cell2/naked2','Cell1/cell2')


subplot(2,2,4)
hold on
for plot_index = [1 2 3 6]
    whole_ratio = [];
    for chromosome_index = 1:length(yeast_genome)
        whole_ratio = [whole_ratio mutation_ratio{plot_index}.penta_norm{chromosome_index}];
    end
    
    ratio_sorted = sort(whole_ratio,'ascend','MissingPlacement','last');
    plot(ratio_sorted,'LineWidth',2);
end
yline(0,'LineWidth',1)
xlim([length(find(~isnan(ratio_sorted)))-10000 length(find(~isnan(ratio_sorted)))])
ylim([-7 7])
ylabel('Cell/naked deam-seq ratio (log2)')
legend('Cell/naked pooled','Cell1/naked1','Cell2/naked2','Cell1/cell2')

