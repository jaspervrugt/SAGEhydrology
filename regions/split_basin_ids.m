function outputlist = split_basin_ids(idlist,K_t,K_v,seed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPLIT_BASIN_IDS Randomize and split basin ID list into train/evaluation
%
% SYNOPSIS:
%   outputlist = split_basin_ids(idlist,K_t,K_v,seed)
%
%   idlist      basin-ID list file name or full path
%               Examples:
%                 '499_basins.txt'
%                 '531_basins.txt'
%                 '671_basins.txt'
%               File must contain basin IDs in a row or column vector
%   K_t         # training basins
%   K_v         # evaluation basins
%   seed        random seed for reproducible shuffling
%               Example: 1, 20260319, etc.
%   outputlist  OUTPUT: string column vector with:
%                 - first K_t IDs sorted in increasing order
%                 - last K_v IDs sorted in increasing order
%
% Note: Writes output file to basins directory, alongside the original file
%       split_xxx_Kt_Kv_seed.txt
%       where xxx is taken from the original file name, e.g.
% Example: outputlist = split_basin_ids('671_basins.txt',500,171,1);
%          -> split_671_500_171_1.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Written for Jasper A. Vrugt, Mar. 2026                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basin directory. The output file is written next to the input file.
dir_basins = pwd;

% Checks
if nargin ~= 4
    error('split_basin_ids:nargin', ...
        ['Function requires 4 inputs: ' ...
         'idlist, K_t, K_v, and seed.']);
end

if ~(isscalar(K_t) ...
        && isnumeric(K_t) ...
        && isfinite(K_t) ...
        && K_t >= 0 ...
        && K_t == floor(K_t))
    error(['      Error: split_basin_ids:K_t', ...
        'K_t must be a nonnegative integer scalar.']);
end

if ~(isscalar(K_v) ...
        && isnumeric(K_v) ...
        && isfinite(K_v) ...
        && K_v >= 0 ...
        && K_v == floor(K_v))
    error(['      Error: split_basin_ids:K_v', ...
        'K_v must be a nonnegative integer scalar.']);
end

if ~(isscalar(seed) ...
        && isnumeric(seed) ...
        && isfinite(seed) ...
        && seed == floor(seed))
    error(['      Error: split_basin_ids:seed', ...
        'seed must be an integer scalar.']);
end

% Resolve input file
if isstring(idlist) || ischar(idlist)
    idlist = char(idlist);
else
    error(['      Error: split_basin_ids:idlist', ...
        'idlist must be a file name or full file path.']);
end

if exist(idlist,'file') == 2
    file_in = idlist;
else
    file_in = fullfile(dir_basins,idlist);
end

if exist(file_in,'file') ~= 2
    error(['      Error: split_basin_ids:fileNotFound', ...
        'Could not find basin list file:\n%s'],file_in);
end

% Read basin IDs. Keep alphanumeric basin IDs intact.
%
% Older versions extracted only digit tokens using regexp(txt,'\d+'),
% which changed IDs such as DE110000, DEA10000 and DEB10000 to numeric
% values and then printed them as 00110000, 00100000, etc. CAMELS-DE,
% CAMELS-AU and other regional datasets require the character prefix to be
% preserved. Purely numeric IDs are still padded to 8 digits to keep the
% original CAMELS-US behavior.
txt = fileread(file_in);

tok = regexp(txt,'[^\s,;]+','match');
tok = string(tok(:));
tok = strip(tok);
tok(tok == "") = [];

if isempty(tok)
    error(['      Error: split_basin_ids:noIDs', ...
        'No basin IDs found in file:\n%s'],file_in);
end

ids = strings(numel(tok),1);
for i = 1:numel(tok)
    s = upper(strip(tok(i)));
    s = erase(s,'"');
    s = erase(s,'''');

    % For pure numeric IDs, retain the existing 8-digit CAMELS-US format.
    if ~isempty(regexp(char(s),'^\d+$','once'))
        ids(i) = string(sprintf('%08d',str2double(s)));
    else
        % For DE110000, DEA10000, AU IDs such as 912101A, etc., preserve
        % the complete alphanumeric identifier.
        ids(i) = s;
    end
end

% Check uniqueness and size
n = numel(ids);

if numel(unique(ids)) ~= n
    warning(['      Warning: split_basin_ids:duplicates', ...
        'Input list contains duplicate basin IDs.']);
end

if K_t + K_v ~= n
    error(['      Error: split_basin_ids:badSplit', ...
        'K_t + K_v = %d, but file contains %d basin IDs.'], ...
        K_t + K_v,n);
end

% Randomize with user-provided seed
rng(seed,'twister');

p = randperm(n);
ids_rand = ids(p);

% Split and sort each subgroup
group_t = sort(ids_rand(1:K_t));
group_v = sort(ids_rand(K_t+1:K_t+K_v));

outputlist = [group_t; group_v];

% Build output file name
[in_dir,name,~] = fileparts(file_in);
if isempty(in_dir)
    in_dir = dir_basins;
end

% Extract leading number from file name, e.g. 499 from 499_basins
tokname = regexprep(name,'_basins$','');

% For names such as CL_297_basins, keep CL_297.
% For names such as 671_basins, keep 671.
% For names without _basins, keep the full name.

file_out = fullfile(in_dir, ...
    sprintf('split_%s_%d_%d_%d.txt', ...
    tokname,K_t,K_v,seed));

% Write output file
writelines(outputlist,file_out);

% Print summary
fprintf('\n');
fprintf('split_basin_ids summary\n');
fprintf('  Input file   : %s\n',file_in);
fprintf('  Total basins : %d\n',n);
fprintf('  K_t          : %d\n',K_t);
fprintf('  K_v          : %d\n',K_v);
fprintf('  Seed         : %d\n',seed);
fprintf('  Output file  : %s\n',file_out);
fprintf('\n');

end