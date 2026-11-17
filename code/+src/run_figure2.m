function run_figure2(out_dir)
%RUN_FIGURE2 Recompute the first-order structure test from initial data.
if nargin<1,error('Specify an output directory for newly computed data.');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
par=model.structure_parameters;
init=model.structure_initial;
sol=src.ms_tumor1d_solver(par,init);
x=sol.x;history=sol.history;diagnostics=sol.diagnostics;parameters=par;
save(fullfile(out_dir,'structure_history.mat'),'x','history','diagnostics','parameters','-v7');
fprintf('Finished the fixed-step Figure 2 computation.\n');
end
