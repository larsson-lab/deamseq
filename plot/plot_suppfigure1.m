%% Load data
load ../workspace.mat;
load ../extra_data.mat;


%% Calculate mean damage levels (allele frequencies) in 250 bp genomic windows
%  (Binned by DNase accessibility)
w_dnase = [];
w_damage_cell = [];
w_damage_naked = [];
for chr_idx = 1:length(yeast_genome)
    idx = strcmp(dnase_tagcounts.chr,yeast_genome(chr_idx).Header);
    
    clear dnase;
    dnase(dnase_tagcounts.start(idx) + 1) = dnase_tagcounts.cleavage(idx); % Positions are zero-based on the DNase dataset, hence +1
    
    f_cell = f{chr_idx}(10,:);
    f_naked = f{chr_idx}(11,:);
    idx_ok = idx_cov_ok{chr_idx};
    
    wsize = 250;
    for i = 1:(floor(min(length(f_cell),length(dnase))/wsize)-1)
        wstart = ((i-1)*wsize+1);
        wend = i*wsize;
        if sum(idx_ok(wstart:wend)) == wsize % Only consider windows that are completely whitelisted
            w_dnase(end+1) = sum(dnase(wstart:wend),'omitnan');
            w_damage_cell(end+1) = mean(f_cell(wstart:wend),'omitnan');
            w_damage_naked(end+1) = mean(f_naked(wstart:wend),'omitnan');
        end
    end
end


%% Partition windows into bins of approx. equal size based on DNase counts
n_bins = 10;

l = quantile(w_dnase,n_bins - 1);
n = discretize(w_dnase,[-inf l inf]);

subplot(3,1,1);
boxplot(w_damage_cell*100,n);
axis([.5 n_bins+0.5 0 1]);
title('Deam-seq, cellular DNA');
ylabel('Mean damage (% allele frequency)');
xlabel('DNase accessibility (low to high)');

subplot(3,1,2);
boxplot(w_damage_naked*100,n);
axis([.5 n_bins+0.5 0 1]);
title('Deam-seq, naked DNA');
ylabel('Mean damage (% allele frequency)');
xlabel('DNase accessibility (low to high)');

subplot(3,1,3);
boxplot(w_dnase,n);
axis([.5 n_bins+0.5 0 1500]);
title('DNase-seq');
ylabel('Total tagcount');
xlabel('DNase accessibility (low to high)');


