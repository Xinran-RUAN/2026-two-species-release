function paths=startup_paper(data_root)
%STARTUP_PAPER Make the three MATLAB packages available from one code path.
% PATHS = STARTUP_PAPER returns the code, release and data locations.
% STARTUP_PAPER(DATA_ROOT) selects a different numerical data tree.
if nargin<1,data_root=[];end
code_root=fileparts(mfilename('fullpath'));addpath(code_root);
paths=utils.paths(data_root);
end
