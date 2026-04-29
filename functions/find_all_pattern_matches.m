function start_index = find_all_pattern_matches(str,pattern)
% Finds the start index of all matches for pattern in str using regexp
% twice to also handle overlapping matches.

start_index = regexp(str,pattern);

% Find overlapping matches that regexp misses:
altered_str = str;
altered_str(start_index) = 'N';
start_index = [start_index regexp(altered_str,pattern)];
