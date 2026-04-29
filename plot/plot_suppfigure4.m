addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat


%% Upstream...
upstream_genes = genes;
upstream_genes.start(upstream_genes.strand == '+') = upstream_genes.start(upstream_genes.strand == '+') - 200;
upstream_genes.start(upstream_genes.strand == '-') = upstream_genes.start(upstream_genes.strand == '-') + 200;
upstream_genes.stop = upstream_genes.start;

ratio_index = [6 2 3];

% Upstream:
ratio_upstream = cell(1,length(ratio_index));
for index = 1:length(ratio_index)
    ratio_upstream{index} = compute_property_around_sites(yeast_genome,mutation_ratio{ratio_index(index)}.penta_norm,200,upstream_genes);
end


%% Downstream...
downstream_genes = genes;
downstream_genes.start(downstream_genes.strand == '+') = downstream_genes.start(downstream_genes.strand == '+') + 200;
downstream_genes.start(downstream_genes.strand == '-') = downstream_genes.start(downstream_genes.strand == '-') - 200;
downstream_genes.stop = downstream_genes.start;

% Downstream:
ratio_downstream =cell(1,length(ratio_index));
for index = 1:length(ratio_index)
    ratio_downstream{index} = compute_property_around_sites(yeast_genome,mutation_ratio{ratio_index(index)}.penta_norm,200,downstream_genes);
end


%% Find max. abs signal of each gene
% Upstream:
strongest_signal = cell(1,length(ratio_index));
for index = 1:length(ratio_index)
    strongest_signal{index} = nan(1,length(upstream_genes.name));
    for gene_index = 1:length(upstream_genes.name)
        idx = abs(ratio_upstream{index}(:,gene_index)) == max(abs(ratio_upstream{index}(:,gene_index)));
        if ~isempty(find(idx,1))
            strongest_signal{index}(gene_index) = ratio_upstream{index}(idx,gene_index);
        end
    end
end

figure()
hold on
for plot_index = 2:length(ratio_index)
    plot(sort(abs(strongest_signal{1})),sort(abs(strongest_signal{plot_index})),'.')
end
line([0 8],[0 8],'Color','k')
axis([0 7 0 7])
xlabel('cell1/cell2')
legend('cell1/naked1','cell2/naked2')
title('Upstream')


% Downstream:
strongest_signal = cell(1,length(ratio_index));
for index = 1:length(ratio_index)
    strongest_signal{index} = nan(1,length(downstream_genes.name));
    for gene_index = 1:length(downstream_genes.name)
        idx = abs(ratio_downstream{index}(:,gene_index)) == max(abs(ratio_downstream{index}(:,gene_index)));
        if ~isempty(find(idx,1))
            strongest_signal{index}(gene_index) = ratio_downstream{index}(idx,gene_index);
        end
    end
end

figure()
hold on
for plot_index = 2:length(ratio_index)
    plot(sort(abs(strongest_signal{1})),sort(abs(strongest_signal{plot_index})),'.')
end
line([0 8],[0 8],'Color','k')
axis([0 7 0 7])
xlabel('cell1/cell2')
legend('cell1/naked1','cell2/naked2')
title('Downstream')


%% All genomic positions:
ratio_all = cell(1,length(ratio_index));
for index = 1:length(ratio_index)
    ratio_all{index} = [];
    for chromosome_index = 1:length(yeast_genome)
        ratio_all{index} = [ratio_all{index} mutation_ratio{ratio_index(index)}.penta_norm{chromosome_index}];
    end
end

figure()
hold on
for plot_index = 2:length(ratio_index)
    plot(sort(abs(ratio_all{1})),sort(abs(ratio_all{plot_index})),'.')
end
line([0 8],[0 8],'Color','k')
axis([0 7 0 7])
xlabel('cell1/cell2')
legend('cell1/naked1','cell2/naked2')
title('Whole-genome')

