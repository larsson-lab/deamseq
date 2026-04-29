addpath('./functions');

%% Load reference genome for Saccaromyces cerevisiae 3:
yeast_genome = fastaread('data/sacCer3/sacCer3.fa');
yeast_genome = yeast_genome(~strcmp({yeast_genome.Header},'chrM')); % Disregard chrM (low mappability)


%% Load and process Deam-seq data
files = {
    'data/calls/IM16/IM16-CELLS-1_dedup' % IM16 = UVC on ice
    'data/calls/IM16/IM16-CELLS-2_dedup'
    'data/calls/IM16/IM16-NAKED-1_dedup'
    'data/calls/IM16/IM16-NAKED-2_dedup'
    'data/calls/IM16/IM16-AA-CELLS_dedup'
    'data/calls/IM16/IM16-AA-CELLS-RAP_dedup'
    'data/calls/IM16/IM16-AA-NAKED_dedup'
    'data/calls/IM8/NO_UV_dedup'        % No UV
    'data/calls/IM9/NO_UV_dedup'        
    }; % The name of all (.tsv)-files with mutation calls

name = {
    'Cells 1'
    'Cells 2'
    'Naked 1'
    'Naked 2'
    'AA Cells'
    'AA Cells + Rap'
    'AA Naked'
    'No UV 1'
    'No UV 2'
    'Cells pooled' % pooled data to be computed further below
    'Naked pooled'
    'No UV pooled'
    }; % Sample names corresponding to the files

% Load the mutation calls (to format cell{chr}(file,pos))
mut = load_mutation_calls(yeast_genome,files);

% Annotate each sample & position with a substitution type
mut = add_substitution_type(yeast_genome,mut);


%% Pool replicate samples (appended to bottom of mut)

for chr_idx = 1:length(yeast_genome)
    mut.n{chr_idx}(length(files)+1:length(files)+3,:) = mut.n{chr_idx}([1 3 8],:) + mut.n{chr_idx}([2 4 9],:);
    mut.cov{chr_idx}(length(files)+1:length(files)+3,:) = mut.cov{chr_idx}([1 3 8],:) + mut.cov{chr_idx}([2 4 9],:);
end


%% Coverage stats

for sample_index = 1:length(name)
    cov = [];
    for j = 1:length(mut.cov)
        cov = [cov double(mut.cov{j}(sample_index,:))];
    end
    mean_cov(sample_index) = mean(cov);
    total_bp(sample_index) = sum(cov);
end

clear cov;

%% Define whitelisted positions with good coverage in all samples
% Exclude overage >=5000 and <=50000 in all main UV samples, and >=3000 in the
% pooled no UV control. The no pooled UV control has uneven coverage -
% hence no upper coverage limit

for chr_idx = 1:length(yeast_genome)
    idx_cov_ok{chr_idx} = sum(mut.cov{chr_idx}(1:4,:) >= 5000 & mut.cov{chr_idx}(1:4,:) <= 50000) == 4 & mut.cov{chr_idx}(12,:) >= 3000;
    idx_cov_ok{chr_idx} = (smooth(idx_cov_ok{chr_idx},500) > 0.9999)'; % expand/filter the inclusion flag
end

%% Calculate mutation fractions (VAFs)
% Calculate per-sample mutation fractions, retaining non-c-diPy positions and preexisting variants

f_all = compute_mutation_fraction(mut);

% Exclude non-whitelisted regions
for chr_idx = 1:length(yeast_genome)
    f_all{chr_idx}(:,~idx_cov_ok{chr_idx}) = NaN;
end


%% Exclude non-diPy positions
idx_c_dipy = find_PyCs(yeast_genome);

f = f_all;
for idx = 1:length(yeast_genome)
    f{idx}(:,~idx_c_dipy{idx}) = NaN;
end


%% Remove preexisting variants
% Defining prexisting variants as >1% VAF present in the pooled No UV
% sample (the flag also includes non-PyC positions)

% non-AA samples
for chr_idx = 1:length(yeast_genome)
    idx_preexisting{chr_idx} = f_all{chr_idx}(12,:) > 0.01;
    f{chr_idx}([1:4 8:12],idx_preexisting{chr_idx}) = NaN;
end

% AA samples (different strain)
for chr_idx = 1:length(yeast_genome)
    idx_preexisting_AA{chr_idx} = f_all{chr_idx}(7,:) > 0.05; % >5% in the naked AA sample
    f{chr_idx}(5:7, idx_preexisting_AA{chr_idx}) = NaN;
end


%% Find cytosine-centered, CPD-compatible pentanucleotides in the genome
[penta_idx,penta_CPD_compatible] = find_CPDpenta_indices(yeast_genome);


%% Calculate sample/sample ratios, normalized by pentanucleotide sequence contexts
% Suggestion: add pooled data last in mut above, to avoid using a separate
% struct for pooled data.

ratio_names = {
    'Cells vs naked pooled'
    'Cells vs naked 1'
    'Cells vs naked 2'
    'AA cells vs naked'
    'AA cells+Rap vs naked'
    'Cells 1 vs cells 2'
    };

mutation_ratio{1}.unnorm = compute_mutation_ratio(f,10,11);
mutation_ratio{2}.unnorm = compute_mutation_ratio(f,1,3);
mutation_ratio{3}.unnorm = compute_mutation_ratio(f,2,4);
mutation_ratio{4}.unnorm = compute_mutation_ratio(f,5,7);
mutation_ratio{5}.unnorm = compute_mutation_ratio(f,6,7);
mutation_ratio{6}.unnorm = compute_mutation_ratio(f,1,2);

for ratio_idx = 1:length(ratio_names)
    % Normalize mutation frequency ratio by subtracting pentanucleotide medians:
    mutation_ratio{ratio_idx}.penta_norm = normalize_ratio_by_penta_median(mutation_ratio{ratio_idx}.unnorm,penta_idx);
    
    % Normalize mutation frequency ratio by subtracting total median:
    mutation_ratio{ratio_idx}.depth_norm = normalize_ratio_by_median(mutation_ratio{ratio_idx}.unnorm);
end


%% Write wig files
write_wig('wig/Cells_naked_pooled',mutation_ratio{1}.penta_norm,yeast_genome);
write_wig('wig/Cells_naked_1',mutation_ratio{2}.penta_norm,yeast_genome);
write_wig('wig/Cells_naked_2',mutation_ratio{3}.penta_norm,yeast_genome);
write_wig('wig/AA_cells_vs_naked',mutation_ratio{4}.penta_norm,yeast_genome);
write_wig('wig/AA_cells_RAP_vs_naked',mutation_ratio{5}.penta_norm,yeast_genome);


%% Cache the workspace
clear ratio_idx chr_idx idx cov j;
save('workspace',"-v7.3");

