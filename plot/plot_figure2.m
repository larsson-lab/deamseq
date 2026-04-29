%% Load data

load ../workspace.mat
addpath('../functions');

set(0, 'DefaultFigureRenderer', 'painters');
figure('Position',[200 500 1000 200]);

palette_16 = readmatrix('../colormaps/palette_16','FileType','text');
palette_22 = readmatrix('../colormaps/palette_22','FileType','text');


%% stacked bar graph with substitution types

sample_order = [1 2 3 4 8 9];

% diPy C positions:
ax = subplot(1,5,1);
cla;

n_mut = zeros(length(files),7);
for chr_idx = 1:length(yeast_genome)
    for sample_idx = 1:length(files)
        for subst_idx = 1:7
            idx = (mut.subst{chr_idx}(sample_idx,:) == subst_idx) & idx_cov_ok{chr_idx} & idx_c_dipy{chr_idx} & ~idx_preexisting{chr_idx}; % pre-existing or not make no practical difference in the plot - for simplicity, excluded from the start
            n_mut(sample_idx,subst_idx) = n_mut(sample_idx,subst_idx) + sum(mut.n{chr_idx}(sample_idx,idx));
        end
    end
end

b = bar(n_mut(sample_order,end:-1:1),'stacked','linestyle','none'); % subst types in reverse order; group samples by experimental condition
lg_labels = {'C>A' 'C>G' 'C>T' 'T>A' 'T>C' 'T>G' 'Indel'};
legend(lg_labels(end:-1:1));
axis([0 7 0 2e8]);
ax.XAxis.TickLength = [0 0];
xticklabels(name(sample_order))
ylabel('Total substitutions');
title('PyC sites')
axis square;
colors = palette_16([6 15 14 11 10 9 7],:);
for idx = 1:length(b)
    b(idx).FaceColor = colors(idx,:);
end


% non-diPy C positions:
ax = subplot(1,5,2);
cla;

n_mut = zeros(length(files),7);
for chr_idx = 1:length(yeast_genome)
    for sample_idx = 1:length(files)
        for subst_idx = 1:7
            idx = (mut.subst{chr_idx}(sample_idx,:) == subst_idx) & idx_cov_ok{chr_idx} & ~idx_c_dipy{chr_idx} & ~idx_preexisting{chr_idx};
            n_mut(sample_idx,subst_idx) = n_mut(sample_idx,subst_idx) + sum(mut.n{chr_idx}(sample_idx,idx));
        end
    end
end

b = bar(n_mut(sample_order, end:-1:1), 'stacked', 'linestyle', 'none'); % subst types in reverse order; group samples by experimental condition
lg_labels = {'C>A' 'C>G' 'C>T' 'T>A' 'T>C' 'T>G' 'Indel'};
legend(lg_labels(end:-1:1));
axis([0 7 0 2e8]);
ax.XAxis.TickLength = [0 0];
xticklabels(name(sample_order))
ylabel('Total substitutions');
title('Other positions')
axis square;
colors = palette_16([6 15 14 11 10 9 7],:);
for idx = 1:length(b)
    b(idx).FaceColor = colors(idx,:);
end


% Number of diPy-C and other positions in whitelisted regions
n1 = 0;
n2 = 0;
for i = 1:16
    n1 = n1 + sum(idx_c_dipy{i} & idx_cov_ok{i}  & ~idx_preexisting{i});
    n2 = n2 + sum(~idx_c_dipy{i} & idx_cov_ok{i} & ~idx_preexisting{i});
end
disp(n1)
disp(n2)


%% Cumulative assayable genome (C-diPy) with frequency above x
% Pooled samples
n_genome = cell(1,length(name));
f_genome = cell(1,length(name));
for sample_idx = 1:length(name) 
    n_genome{sample_idx} = [];
    f_genome{sample_idx} = [];
    for chr_idx = 1:length(yeast_genome)
        n_genome{sample_idx} = [n_genome{sample_idx} mut.n{chr_idx}(sample_idx,idx_c_dipy{chr_idx} & idx_cov_ok{chr_idx} & ~idx_preexisting{chr_idx})]; % preexisting makes little difference
        f_genome{sample_idx} = [f_genome{sample_idx} f{chr_idx}(sample_idx,:)]; % preexisting, non-PyC, or low-coverage positions already gone
    end
end

n_range = 0:200;
sample_order = [10 11 12];
disp(name(sample_order));

clear gen_frac_n;
for j = 1:length(sample_order)
    for i = 1:length(n_range)
        gen_frac_n(j, i) = mean(n_genome{sample_order(j)} >= n_range(i));
    end
end

ax = subplot(1,5,3);
cla

colororder(ax, palette_22([9 7 5 1],:));
hold on;

for idx = 1:length(sample_order)
    plot(n_range,gen_frac_n(idx,:),'LineWidth',1,'Color',palette_22(idx,:));
end
ylabel('Fraction of diPy-C pos.');
xlabel('Substitutions (n >=)');
legend({'Cells' 'Naked' 'No UV'});
axis square;


disp(gen_frac_n(:,n_range == 20)) % report for paper
disp(gen_frac_n(:,n_range == 60)) % report for paper

% VAF/n stats
disp(mean(f_genome{10}*100,'omitnan')); % cells pooled
disp(mean(f_genome{11}*100,'omitnan')); % naked pooled
disp(mean(f_genome{12}*100,'omitnan')); % no UV pooled

disp(mean(n_genome{10},'omitnan'));
disp(mean(n_genome{11},'omitnan'));
disp(mean(n_genome{12},'omitnan'));


%% Scatterplot - downsampled to 0.5%

ax = subplot(1,5,5);
cla
hold on;

dstep = 200;

f_genome = cell(1,length(name));
for sample_idx = 1:length(name)
    f_genome{sample_idx} = [];
    for chr_idx = 1:length(yeast_genome)
        f_genome{sample_idx} = [f_genome{sample_idx} f{chr_idx}(sample_idx, :)]; 
    end
end

scatter(f_genome{1}(1:dstep:end)*100, f_genome{12}(1:dstep:end)*100,75,'.'); % Downsample to 0.5%
scatter(f_genome{1}(1:dstep:end)*100, f_genome{11}(1:dstep:end)*100,75,'.');
scatter(f_genome{1}(1:dstep:end)*100, f_genome{2}(1:dstep:end)*100,75,'.');
axis([0 3 0 3]);
axis square;
xticks([0 1 2 3]);
yticks([0 1 2 3]);
xlabel('Cells 1')
ylabel('Sync. 2 Unsync. 1 Naked No UV')
colororder(ax, palette_22([1 5 7 9],:))
legend({'No UV vs cells 1' 'Naked vs cells 1' 'Cells 2 vs cells 1'});

% Pearson correlation:
corr(f_genome{1}(~isnan(f_genome{1}))', f_genome{12}(~isnan(f_genome{4}))')
corr(f_genome{1}(~isnan(f_genome{1}))', f_genome{11}(~isnan(f_genome{3}))')
corr(f_genome{1}(~isnan(f_genome{1}))', f_genome{2}(~isnan(f_genome{2}))')
name([12 11 2])


%% Prepare for stacked VAF trinucleotides histogram plot

[tri,tri_rv] = find_trinucleotides(yeast_genome);


%% Stacked VAF histogram - per trinucleotide
% preexisting variants not considered (based on f)

trinuc = {'CCC','CCG','CCA','CCT','TCC','TCT','TCA','TCG','ACC','GCC','ACT','GCT'};
trinuc = trinuc([10 12 2 8 7 9 11 3 1 4 5 6]);

hist_range = 0:0.0001:0.02;
i_sample = 10; % pooled cells

ax = subplot(1,5,4);
cla;

clear stacked_hist;
for tri_idx = 1:length(trinuc)
    disp(tri_idx)
    stacked_hist(tri_idx,:) = zeros(size(hist_range));
    for chr_idx = 1:length(yeast_genome)
        idx = strcmp(tri{chr_idx},trinuc{tri_idx}) | strcmp(tri_rv{chr_idx},trinuc{tri_idx});
        x = f{chr_idx}(i_sample,idx);
        x(isnan(x) | x > max(hist_range)) = [];
        stacked_hist(tri_idx,:) = stacked_hist(tri_idx,:) + hist(x,hist_range);
    end
end

area(hist_range*100,stacked_hist','LineStyle','none');
xlabel('Variant allele freq. (%)');
ylabel('DiPy-C pos. (n)');
legend(trinuc);


% superimpose the non-exposed pooled sample - using only one color
i_sample = 12;

nonstacked_hist = zeros(size(hist_range));
for chr_idx = 1:length(yeast_genome)
    x = f{chr_idx}(i_sample,:);
    x(isnan(x) | x > max(hist_range)) = [];
    nonstacked_hist = nonstacked_hist + hist(x,hist_range);
end
hold on;
area(hist_range*100,nonstacked_hist,'LineStyle','none','FaceColor',[.6 .6 .6]);
axis([0 2 0 7e4]);
axis square;
colororder(ax,palette_22(5:16,:));


%% example region
% NOTE: all positions are shown including non-PyC, which illustrates good signal to
% noise
% NOTE2: no prexisting variants in this region

figure('Position',[200 500 1000 380]);
colororder(palette_22(3:end,:));

chr = 1;
region = 50000:50120;
idx_cov_ok{chr}(region)

sample_order = [1 3 8];

for i = 1:length(sample_order)
    subplot(5,1,i);
    bar(region,[mut.n{chr}(sample_order(i),region); mut.n{chr}(sample_order(i)+1,region)], 'stacked', 'LineStyle','none');
    xticks([]);
    text(region(1),400,name{sample_order(i)});
    axis([region(1)-1 region(end)+1 0 400]);
end

for idx = 1:length(region)
    if idx_c_dipy{chr}(region(idx))
        c = 'k';
    else
        c = 'y';
    end
    text(region(idx), -100, yeast_genome(chr).Sequence(region(idx)),'HorizontalAlignment','center', 'Color', c);
    text(region(idx), -200, seqrcomplement(yeast_genome(chr).Sequence(region(idx))),'HorizontalAlignment','center', 'Color', c);
end

text(region((end+1)/2), 2650, [yeast_genome(chr).Header ':' num2str(region(1)) '-' num2str(region(end))],'HorizontalAlignment','center');


%% report some numbers for paper

% mutation frequency
f_genome = [];
for chr_idx = 1:length(yeast_genome)
    f_genome = [f_genome f{chr_idx}];
end
disp(mean(f_genome, 2, 'omitnan'))

% mutation density (mutations per bp sequenced), pooled samples, same order
for sample_idx = 1:length(name)
    n_mutations = 0;
    n_sequenced = 0;
    
    for i = 1:length(yeast_genome)
        n_mutations = n_mutations + sum(mut.n{i}(sample_idx,idx_c_dipy{i} & idx_cov_ok{i} & ~idx_preexisting{i})); % total mutations mapped, disregarding those at non-diPy postions, only in included regions unlikely to have mapping artefacts
        n_sequenced = n_sequenced + sum(mut.cov{i}(sample_idx,idx_cov_ok{i} & ~idx_preexisting{i})); % total bp sequenced, also only in included regions
    end

    fprintf('%s\t%1.2f\n', name{sample_idx}, n_sequenced/n_mutations);
end

