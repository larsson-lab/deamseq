addpath('../functions');

%% Load data
load ../workspace.mat
load ../extra_data.mat

% Load gene annotation data:
annot_gtf = GTFAnnotation('../data/sacCer3/sacCer3.ncbiRefSeq.gtf');
annot_genes = annot_gtf.getGenes;
annot_genes = annot_genes(~strcmp(string(annot_genes.Reference),'chrM'),:);

% Load reference genome for Saccaromyces cerevisiae 3:
yeast_genome = fastaread('../data/sacCer3/sacCer3.fa');
yeast_genome = yeast_genome(~strcmp({yeast_genome.Header},'chrM')); % Disregard chrM (low mappability)

% Load motif sequences for transcription factors:
motifs = readtable('../data/motifs/TF_motif_sequences_short.txt','FileType','text','ReadVariableNames',false,'Delimiter',',');
motifs = renamevars(motifs,["Var1","Var2","Var3","Var4"],["Name","Seq","Sym_seq","Sym_seq_rcomp"]);


%% Load RNA-seq data for Reb1 depletion
expr_plus = read_wig('../data/REB1AA_RNAseq/GSM2580262_Reb1-AARrp6delta-Rapa_plus.wig');
expr_minus = read_wig('../data/REB1AA_RNAseq/GSM2580262_Reb1-AARrp6delta-Rapa_minus.wig');
expr_plus_rap = read_wig('../data/REB1AA_RNAseq/GSM2580263_Reb1-AARrp6delta+Rapa_plus.wig');
expr_minus_rap = read_wig('../data/REB1AA_RNAseq/GSM2580263_Reb1-AARrp6delta+Rapa_minus.wig');


expr_plus(strcmp(expr_plus.Chromosome,'2-micron'),:) = [];
expr_minus(strcmp(expr_minus.Chromosome,'2-micron'),:) = [];
expr_plus_rap(strcmp(expr_plus_rap.Chromosome,'2-micron'),:) = [];
expr_minus_rap(strcmp(expr_minus_rap.Chromosome,'2-micron'),:) = [];

expr_plus(strcmp(expr_plus.Chromosome,'chrMito'),:) = [];
expr_minus(strcmp(expr_minus.Chromosome,'chrMito'),:) = [];
expr_plus_rap(strcmp(expr_plus_rap.Chromosome,'chrMito'),:) = [];
expr_minus_rap(strcmp(expr_minus_rap.Chromosome,'chrMito'),:) = [];


%% Compute per-gene expression as counts per gene
% -RAP:
expression = nan(1,length(annot_genes.GeneID));
for gene_idx = 1:length(annot_genes.GeneID)
    gene_region = annot_genes.Start(gene_idx):annot_genes.Stop(gene_idx);

    if annot_genes.Strand(gene_idx) == '+'
        current_chr = strcmp(string(annot_genes.Reference(gene_idx)),expr_plus.Chromosome);

        expr_plus_chr = expr_plus(current_chr,:);
        expr_region = expr_plus_chr.Value(ismember(expr_plus_chr.Position,gene_region));

    elseif annot_genes.Strand(gene_idx) == '-'
        current_chr = strcmp(string(annot_genes.Reference(gene_idx)),expr_minus.Chromosome);

        expr_minus_chr = expr_minus(current_chr,:);
        expr_region = expr_minus_chr.Value(ismember(expr_minus_chr.Position,gene_region));
    end

    expression(gene_idx) = sum(expr_region);
end


% +RAP:
expression_rap = nan(1,length(annot_genes.GeneID));
for gene_idx = 1:length(annot_genes.GeneID)
    gene_region = annot_genes.Start(gene_idx):annot_genes.Stop(gene_idx);

    if annot_genes.Strand(gene_idx) == '+'
        current_chr = strcmp(string(annot_genes.Reference(gene_idx)),expr_plus_rap.Chromosome);
        expr_plus_chr = expr_plus_rap(current_chr,:);
        expr_region = expr_plus_chr.Value(ismember(expr_plus_chr.Position,gene_region)); 
    elseif annot_genes.Strand(gene_idx) == '-'
        current_chr = strcmp(string(annot_genes.Reference(gene_idx)),expr_minus_rap.Chromosome);
        expr_minus_chr = expr_minus_rap(current_chr,:);
        expr_region = expr_minus_chr.Value(ismember(expr_minus_chr.Position,gene_region));
    end

    expression_rap(gene_idx) = sum(expr_region);
end

clear expr_plus expr_minus expr_plus_rap expr_minus_rap


%% Compare -RAP and +RAP
expr.name = string(annot_genes.GeneName');
expr.expression = expression;
expr.expression_rap = expression_rap;
expr.ratio = (expr.expression_rap + 1) ./ (expr.expression + 1); % One pseudocount

[~,sort_idx] = sort(expr.ratio,'ascend','MissingPlacement','last');

expr.name = expr.name(sort_idx);
expr.expression = expr.expression(sort_idx);
expr.expression_rap = expr.expression_rap(sort_idx);
expr.ratio = expr.ratio(sort_idx);


%% Look at noise
figure()
plot(log2(expr.expression+expr.expression_rap),log2(expr.ratio),'.')


%% Filter out too lowly expressed genes
remove_idx = expr.expression <= 1000;

fields = fieldnames(expr);
for i = 1:length(fieldnames(expr))
    expr.(fields{i}) = expr.(fields{i})(~remove_idx);
end


%% Find all Reb1 sites
TF_name = 'Reb1';
TF_motif_index = find(strcmp(motifs.Name,TF_name));
motif_sequence = motifs.Seq{TF_motif_index};

% Find motif positions:
site = find_seq_motif_sites(yeast_genome,motif_sequence,motifs.Sym_seq{TF_motif_index});

fprintf(['Number of motif sites for ' TF_name ': %d\n'],length(site.chr));


%% Remove all sites that are not included in the whitelisted regions:
outside_whitelisted = false(1,length(site.chr));
for site_index = 1:length(site.chr)
    current_chr = strcmp(site.chr(site_index),{yeast_genome.Header});
    outside_whitelisted(site_index) = any(idx_cov_ok{current_chr}(site.start(site_index):site.stop(site_index)) == 0);
end
site.chr = site.chr(~outside_whitelisted);
site.start = site.start(~outside_whitelisted);
site.stop = site.stop(~outside_whitelisted);
site.strand = site.strand(~outside_whitelisted);
site.motif = site.motif(~outside_whitelisted);


%% Calculate ratio diff around each predicted TF site
site.site_ratio_aa = compute_property_around_sites(yeast_genome,mutation_ratio{4}.penta_norm,0,site);
site.site_ratio_aa_rap = compute_property_around_sites(yeast_genome,mutation_ratio{5}.penta_norm,0,site);

site.ratio_diff = mean(abs(site.site_ratio_aa),'omitnan') - mean(abs(site.site_ratio_aa_rap),'omitnan');


%% Match genes with the gene names in the expr struct, and find RNA-seq ratio for all genes
[~,idx_expr,idx_genes] = intersect(expr.name,genes.name,'stable');

fields = fieldnames(expr);
for i = 1:length(fieldnames(expr))
    expr.(fields{i}) = expr.(fields{i})(idx_expr);
end

fields = fieldnames(genes);
for i = 1:length(fieldnames(genes))
    genes.(fields{i}) = genes.(fields{i})(idx_genes);
end


%% Check if there are Reb1 sites (bound or not bound) 200 bp upstream of genes
has_reb1 = zeros(1,length(genes.name));
has_reb1_bind = zeros(1,length(genes.name));
has_reb1_bind_weak = zeros(1,length(genes.name));
has_reb1_nobind = zeros(1,length(genes.name));
for gene_idx = 1:length(genes.name)
    sites_current_chr = strcmp(string(genes.chr(gene_idx)),site.chr);

    if genes.strand(gene_idx) == '+'
        upstream_region = genes.start(gene_idx)-200:genes.start(gene_idx);
    elseif genes.strand(gene_idx) == '-'
        upstream_region = genes.start(gene_idx):genes.start(gene_idx)+200;
    end

    has_reb1_temp = ismember(site.start(sites_current_chr),upstream_region) | ismember(site.stop(sites_current_chr),upstream_region);
    has_reb1(gene_idx) = any(has_reb1_temp);

    if has_reb1(gene_idx)
        has_reb1_bind(gene_idx) = any(has_reb1_temp & site.ratio_diff(sites_current_chr) > 0.75);
        has_reb1_bind_weak(gene_idx) = any(has_reb1_temp & site.ratio_diff(sites_current_chr) >= 0.25 & site.ratio_diff(sites_current_chr) < 0.75);
        has_reb1_nobind(gene_idx) = any(has_reb1_temp & site.ratio_diff(sites_current_chr) >= 0 & site.ratio_diff(sites_current_chr) < 0.25);
    end
end


%% Plot
figure()
hold on
cdfplot(log2(expr.ratio(find(has_reb1_bind))));
cdfplot(log2(expr.ratio(find(has_reb1_bind_weak))));
cdfplot(log2(expr.ratio(find(has_reb1_nobind))));
cdfplot(log2(expr.ratio(find(~has_reb1))));
legend(['Strongly bound Reb1 sites (n=' num2str(length(find(has_reb1_bind))) ')'],['Weakly bound Reb1 sites (n=' num2str(length(find(has_reb1_bind_weak))) ')'],['Not bound Reb1 sites (n=' num2str(length(find(has_reb1_nobind))) ')'],['No Reb1 sites (n=' num2str(length(find(~has_reb1))) ')'],'Location','southeast')

