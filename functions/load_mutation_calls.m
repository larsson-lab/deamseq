function mut = load_mutation_calls(ref_genome,files)
% First constructs the struct mut with empty fields n, var and cov based on
% the length of the reference genome .fa-file ref_genome. Then loads the
% tsv files specified by the file names in files, and for
% each file and mutated position sets the values of the mut fields based on
% Reads2 (number of mutated reads), Cov (total coverage) and Var
% (substituted variable) as specified by the tsv files.

for chr_index = 1:length(ref_genome)
    % Note: any diPy position with zero mutations also has low coverage -
    % no need to use NaN or similar to indicate these sites
    mut.n{chr_index} = zeros(length(files),length(ref_genome(chr_index).Sequence),'int32'); % int32 saves 50% of RAM
    mut.cov{chr_index} = zeros(length(files),length(ref_genome(chr_index).Sequence),'int32');
    mut.var{chr_index} = repmat('.',length(files),length(ref_genome(chr_index).Sequence)); % variant base
end

for file_index = 1:length(files)
    disp(file_index)
    calls = readtable([files{file_index} '.tsv'], 'ReadVariableNames', true, 'filetype', 'text', 'Delimiter', '\t');
    idx = find(strcmp(calls.Chrom, 'chrM'));
    calls(idx(1):end, :) = []; % Disregard chrM for this study

    for j = 1:length(calls.Position)
        current_chr = strcmp({ref_genome.Header},calls.Chrom{j});
        current_pos = calls.Position(j);
        if ~strcmp(calls.Ref{j},ref_genome(current_chr).Sequence(current_pos))
           fprintf('Warning! Reference base mismatch at %d (%s)\n',current_pos,calls.Ref{j});
        end
        mut.n{current_chr}(file_index,current_pos) = calls.Reads2(j);
        mut.cov{current_chr}(file_index,current_pos) = calls.Reads1(j) + calls.Reads2(j); % should be a better estimate than the .Cov column
        mut.var{current_chr}(file_index,current_pos) = calls.Var{j}(1); % Indels will be indicated by + or -
    end
end
