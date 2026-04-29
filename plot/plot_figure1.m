addpath('../functions');

load ../workspace.mat
load ../nodedup_stats.mat % contains stats for non-deduplicated samples


%% Show coverage stats
for sample_idx = 1:length(files)
    fprintf('%s\t%1.0f (%1.0f)\t%d (%d)\n',name{sample_idx},mean_cov(sample_idx),mean_cov_nodedup(sample_idx),total_bp(sample_idx),total_bp_nodedup(sample_idx));
end

% Compute the percentage of genome within white-listed regions
included = [];
for chr_idx = 1:length(yeast_genome)
    included = [included idx_cov_ok{chr_idx}];
end
mean(included)
