addpath('./functions');

%% Load data
% Load reference genome for Saccaromyces cerevisiae 3:
yeast_genome = fastaread('data/sacCer3/sacCer3.fa');
yeast_genome = yeast_genome(~strcmp({yeast_genome.Header},'chrM')); % Disregard chrM (low mappability)

% Load gene annotation data:
annot_gtf = GTFAnnotation('data/sacCer3/sacCer3.ncbiRefSeq.gtf');
annot_genes = annot_gtf.getGenes;
annot_genes = annot_genes(~strcmp(string(annot_genes.Reference),'chrM'),:);

% Load gene TSSs from SMORE-seq:
gene_TSS = readtable('data/TSS_data/GSE49026_S-TSS.txt','FileType','text','ReadVariableNames',true,'Delimiter','\t');
gene_TSS = gene_TSS(~strcmp(string(gene_TSS.chr),'chrM'),:);

% Load DNase cleavage and footprints:
% New 2026:
dnase_footprints = readtable('data/DNase/lift_over/sacCer3.yeast.footprints_noheader.bed','FileType','text','ReadVariableNames',false,'Delimiter','tab');
dnase_tagcounts = readtable('data/DNase/lift_over/sacCer3.yeast.dnaseI.tagCounts.bed','FileType','text','ReadVariableNames',false,'Delimiter','tab');
dnase_footprints = renamevars(dnase_footprints,["Var1","Var2","Var3"],["chr","start","stop"]);
dnase_tagcounts = renamevars(dnase_tagcounts,["Var1","Var2","Var3","Var5"],["chr","start","stop","cleavage"]);
dnase_tagcounts = removevars(dnase_tagcounts,{'Var4'});

%% Edit the gene data to only keep the SMORE-seq TSS coordinates
gene_starts = reverse_gene_start_by_strand(annot_genes); % If gene is on reverse strand, change gene start to the gene stop position

genes.start = gene_starts;
edited_genes = [];
for TSS_index = 1:length(gene_TSS.chr)
    index_chr = find(strcmp(string(gene_TSS.chr(TSS_index)),string(annot_genes.Reference)));
    distances_start = abs(gene_starts(index_chr) - gene_TSS.coordinate(TSS_index));

    genes.start(index_chr(distances_start == min(distances_start))) = gene_TSS.coordinate(TSS_index);
    if isscalar(index_chr(distances_start == min(distances_start)))
        edited_genes = [edited_genes index_chr(distances_start == min(distances_start))];
    end
end

genes.name = annot_genes.GeneName(edited_genes);
genes.ID = annot_genes.GeneID(edited_genes);
genes.chr = annot_genes.Reference(edited_genes);
genes.strand = annot_genes.Strand(edited_genes);
genes.start = genes.start(edited_genes);

clear annot_gtf annot_genes gene_TSS gene_starts yeast_genome TSS_index index_chr;


%% Cache the workspace
save('extra_data',"-v7.3");
