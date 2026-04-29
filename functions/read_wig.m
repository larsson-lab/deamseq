function wig_table = read_wig(file_path)
chr_names = {'chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI'}; % TODO

data_table = readtable(file_path,'FileType','text');
variable_step_index = find(isnan(data_table.Var1));

wig_table = table('Size',[0 3],'VariableTypes',{'string','double','double'});

wig_table(1:variable_step_index(1)-1,1) = chr_names(1);
wig_table(1:variable_step_index(1)-1,2:3) = data_table(1:variable_step_index(1)-1,:);
for chr_index = 2:length(chr_names)
    idx = variable_step_index(chr_index-1):variable_step_index(chr_index)-1;

    wig_table(idx,1) = chr_names(chr_index);
    wig_table(idx,2:3) = data_table(idx,:);
end

wig_table(isnan(wig_table.Var2),:) = [];
wig_table = renamevars(wig_table,["Var1","Var2","Var3"],["Chromosome","Position","Value"]);
