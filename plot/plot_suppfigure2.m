addpath('../functions');

%% Load data
load '../workspace.mat'
load ../extra_data.mat

cell_naked_ratio = mutation_ratio{1};

%% Calculate std:
fprintf('Std. for depth norm: %.4f \n',std(cell2mat(cell_naked_ratio.depth_norm),'omitmissing'))
fprintf('Std. for penta norm: %.4f \n',std(cell2mat(cell_naked_ratio.penta_norm),'omitmissing'))


%% Find cytosine-centered, CPD-compatible pentanucleotides in the genome
[~,penta_CPD_compatible] = find_CPDpenta_indices(yeast_genome);


%% Calculate median ratio per pentanucleotide:
median_ratio = compute_penta_median_ratio(cell_naked_ratio.depth_norm,penta_idx);

% Plot median ratio for each pentanucleotide:
figure()
bar(median_ratio);
xticks(1:length(penta_CPD_compatible))
xticklabels(penta_CPD_compatible)
xtickangle(90)
yticks(-0.6:0.2:0.4)
ylim([-0.6 0.4])


%% Plot for a random region:
figure()
a = [cell_naked_ratio.depth_norm{1}(135850:135950); cell_naked_ratio.penta_norm{1}(135850:135950)];
bar(a',1);
xlim([-0.5 100.5])
xticklabels(135850:10:135950)
legend('Depth norm','Penta norm')
