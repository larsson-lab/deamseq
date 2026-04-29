Reb1_chip = readtable('../mmc2.xls','ReadVariableNames',false,'Sheet','Reb1');
Reb1_chip = Reb1_chip(3:end,:);

Reb1_chip = Reb1_chip(:,[1 2 2 3]);


%% (Remove the zeros...)
chr1 = {'chr01','chr02','chr03','chr04','chr05','chr06','chr07','chr08','chr09','chr10','chr11','chr12','chr13','chr14','chr15','chr16'};
chr2 = {'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16'};

new_chr = strings(height(Reb1_chip),1);
for i = 1:height(Reb1_chip)
    new_chr(i) = chr2(strcmp(Reb1_chip.Var1(i),chr1));
end

Reb1_chip.Var1 = new_chr;


%%
writetable(Reb1_chip,'Reb1_peaks.bed','FileType','text','Delimiter','\t','WriteVariableNames',false);
