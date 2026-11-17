function p=paths(data_root)
%PATHS Resolve the code and sibling data tree without using the current folder.
p.code=fileparts(fileparts(mfilename('fullpath')));
p.root=fileparts(p.code);p.data=fullfile(p.root,'data');
if nargin>0 && ~isempty(data_root),p.data=char(data_root);end
p.provenance=fullfile(p.code,'+utils','provenance');
end
