function out_dir = redraw_paper(result, out_dir, data_root)
%REDRAW_PAPER Recreate manuscript figures and tables from pure saved data.
%   REDRAW_PAPER('all', OUT) redraws Figures 1--6 and exports Tables 1--3.
%   REDRAW_PAPER(2, OUT), REDRAW_PAPER('Figure2', OUT), and
%   REDRAW_PAPER('Table1', OUT) select individual results.
%   REDRAW_PAPER('tables', OUT) exports all three tables.
%   REDRAW_PAPER(RESULT, OUT, DATA_ROOT) uses a different data directory.
%   This function never runs a numerical solver and never writes to data/.
if nargin<3,data_root=[];end
paths = startup_paper(data_root);
if nargin < 1 || isempty(result), result = 'all'; end
if nargin < 2 || isempty(out_dir)
    out_dir = fullfile(tempdir, 'two-species-paper-results');
end
if nargin < 3 || isempty(data_root), data_root = paths.data; end
out_dir = char(out_dir);
data_root = char(data_root);
assert(isfolder(data_root), 'two_species:MissingData', ...
    'Data directory does not exist: %s', data_root);
canonical_data = char(java.io.File(data_root).getCanonicalPath());
canonical_out = char(java.io.File(out_dir).getCanonicalPath());
assert(~strcmp(canonical_out,canonical_data) && ...
    ~startsWith(canonical_out,[canonical_data filesep]), ...
    'two_species:ReadOnlyData', 'Choose an output directory outside data/.');
selectors = parse_result(result);
if ~isfolder(out_dir), mkdir(out_dir); end
for k = 1:numel(selectors)
    id = selectors{k};
    assert(isfolder(fullfile(data_root,id)), 'two_species:MissingData', 'Missing result data: %s',id);
    fprintf('Reproducing %s from saved numerical data...\n',id);
    switch id
        case 'Figure1', utils.reproduce_figure1(data_root,out_dir);
        case 'Figure2', utils.reproduce_figure2(data_root,out_dir);
        case 'Figure3', utils.reproduce_figure3(data_root,out_dir);
        case 'Figure4', utils.reproduce_figure4(data_root,out_dir);
        case 'Figure5', utils.reproduce_figure5(data_root,out_dir);
        case 'Figure6', utils.reproduce_figure6(data_root,out_dir);
        case 'Table1', utils.reproduce_table1(data_root,out_dir);
        case 'Table2', utils.reproduce_table2(data_root,out_dir);
        case 'Table3', utils.reproduce_table3(data_root,out_dir);
    end
end
fprintf('Reproduced results are in %s\n', out_dir);
end

function ids = parse_result(result)
if isnumeric(result)
    assert(isvector(result) && ~isempty(result) && ...
        all(isfinite(result)) && all(result == fix(result)) && ...
        all(result >= 1 & result <= 6), 'two_species:InvalidResult', ...
        'Numeric selectors must be figure numbers 1 through 6.');
    ids = arrayfun(@(n) sprintf('Figure%d', n), result, 'UniformOutput', false);
    return
end
assert(ischar(result) || (isstring(result) && isscalar(result)), ...
    'two_species:InvalidResult', 'Use a figure number, FigureN, TableN, or all.');
s = lower(strtrim(char(result)));
if strcmp(s, 'all')
    ids = [arrayfun(@(n) sprintf('Figure%d', n), 1:6, 'UniformOutput', false), ...
        arrayfun(@(n) sprintf('Table%d', n), 1:3, 'UniformOutput', false)];
elseif strcmp(s, 'figures')
    ids = arrayfun(@(n) sprintf('Figure%d', n), 1:6, 'UniformOutput', false);
elseif any(strcmp(s, {'table', 'tables'}))
    ids = arrayfun(@(n) sprintf('Table%d', n), 1:3, 'UniformOutput', false);
elseif ~isempty(regexp(s, '^figure[1-6]$', 'once'))
    ids = {['Figure' s(end)]};
elseif ~isempty(regexp(s, '^table[1-3]$', 'once'))
    ids = {['Table' s(end)]};
else
    error('two_species:InvalidResult', ...
        'Unknown result "%s". Use Figure1--Figure6, Table1--Table3, or all.', s);
end
end
