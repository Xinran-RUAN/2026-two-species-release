function summary = validate_release(data_root)
%VALIDATE_RELEASE Check the saved-data layout and reject executable MAT payloads.
%   This is a read-only packaging check; it does not run a simulation.
if nargin<1,data_root=[];end
paths = startup_paper(data_root);
if nargin < 1 || isempty(data_root), data_root = paths.data; end
expected = [arrayfun(@(n) sprintf('Figure%d',n),1:6,'UniformOutput',false), ...
    arrayfun(@(n) sprintf('Table%d',n),1:3,'UniformOutput',false)];
for k = 1:numel(expected)
    assert(isfolder(fullfile(data_root,expected{k})), ...
        'two_species:MissingData', 'Missing data/%s.', expected{k});
    assert(~isempty(which(['utils.reproduce_' lower(expected{k})])), ...
        'two_species:MissingCode', 'Missing reproduction entry for %s.',expected{k});
end
files = dir(fullfile(data_root, '**', '*'));
files = files(~[files.isdir] & ~startsWith({files.name},'.'));
mat_count = 0; csv_count = 0;
for k = 1:numel(files)
    filename = fullfile(files(k).folder, files(k).name);
    [~,~,extension] = fileparts(filename);
    switch lower(extension)
        case '.mat'
            contents = load(filename);
            check_value(contents,filename);
            mat_count = mat_count + 1;
        case '.csv'
            assert(files(k).bytes > 0, 'two_species:EmptyData', ...
                'Empty CSV: %s',filename);
            csv_count = csv_count + 1;
        case '.json'
            jsondecode(fileread(filename));
        otherwise
            error('two_species:UnexpectedDataFile', ...
                'data/ may contain only numeric MAT, CSV and parameter JSON files: %s', filename);
    end
end
summary = struct('figures',6,'tables',3,'mat_files',mat_count, ...
    'csv_files',csv_count,'total_files',numel(files));
fprintf('Data validation passed: %d figures, %d tables, %d files.\n', ...
    summary.figures,summary.tables,summary.total_files);
end

function check_value(value,location)
if isnumeric(value) || islogical(value) || ischar(value) || isstring(value)
    return
elseif isstruct(value)
    names = fieldnames(value);
    for j = 1:numel(value)
        for k = 1:numel(names)
            check_value(value(j).(names{k}),[location '.' names{k}]);
        end
    end
elseif iscell(value)
    for k = 1:numel(value), check_value(value{k},location); end
else
    error('two_species:ExecutableData', ...
        'Non-numeric payload of class %s in %s.',class(value),location);
end
end
