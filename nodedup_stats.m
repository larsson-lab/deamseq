addpath('./functions');

%% Load sacCer3 genome and define some general stuff
% Load reference genome for Saccaromyces cerevisiae 3:
yeast_genome = fastaread('data/sacCer3/sacCer3.fa');
yeast_genome = yeast_genome(~strcmp({yeast_genome.Header},'chrM')); % Disregard chrM (low mappability)


%% Load and process deamseq CPD data
files = {
    'data/calls/IM16/IM16-CELLS-1'
    'data/calls/IM16/IM16-CELLS-2'
    'data/calls/IM16/IM16-NAKED-1'
    'data/calls/IM16/IM16-NAKED-2'
    'data/calls/IM16/IM16-AA-CELLS'
    'data/calls/IM16/IM16-AA-CELLS-RAP'
    'data/calls/IM16/IM16-AA-NAKED'
    'data/calls/IM8/NO_UV'        % No UV
    'data/calls/IM9/NO_UV'        
    }; % The name of all (.tsv)-files with mutation calls

name = {
    'Cells 1 nodedup'
    'Cells 2 nodedup'
    'Naked 1 nodedup'
    'Naked 2 nodedup'
    'AA Cells nodedup'
    'AA Cells+Rap nodedup'
    'AA Naked nodedup'
    'No UV 1 nodedup'
    'No UV 2 nodedup'    
    }; % Sample names corresponding to the files

% Load the mutation calls (to format cell{chr}(file,pos))
mut = load_mutation_calls(yeast_genome,files);


%% Get coverage stats comparing dedup and non-dedup samples
for sample_index = 1:length(name)
    cov_nodedup = [];
    for j = 1:length(mut.cov)
        cov_nodedup = [cov_nodedup double(mut.cov{j}(sample_index,:))];
    end
    mean_cov_nodedup(sample_index) = mean(cov_nodedup);
    total_bp_nodedup(sample_index) = sum(cov_nodedup);
end


%% Cache workspace
clear files mut name yeast_genome cov_nodedup j;
save ('nodedup_stats');

