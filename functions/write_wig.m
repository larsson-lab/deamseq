function write_wig(file,input_cell,ref_genome)
% Reformats input_cell (cell array where each cell represents one
% chromosome), and writes it as a WIG-file to ../wig/"file_name".wig.

fid = fopen([file '.wig'],'w');
for chr_idx = 1:length(ref_genome)
    header_str = ['variableStep chrom=' ref_genome(chr_idx).Header '\n'];
    fprintf(fid,header_str);

    input_current_chr = input_cell{chr_idx};
    input_current_chr(isnan(input_current_chr)) = 0;
    input_current_chr = string(input_current_chr);

    position = 1:length(input_current_chr);
    position = string(position);

    total_str = strings(1,length(position));
    for pos_idx = 1:length(position)
        total_str(pos_idx) = append(string(position(pos_idx)),'    ',string(input_current_chr(pos_idx)));
        main_str = append(string(position(pos_idx)),'    ',string(input_current_chr(pos_idx)),'\n');
        fprintf(fid,main_str);
    end
end
fclose(fid);
