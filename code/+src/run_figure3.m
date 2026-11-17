function run_figure3(out_dir)
%RUN_FIGURE3 Recompute final-paper two- and three-species structure data.
% The paper default Newton tolerance 1e-12 is used, rather than the archived
% two-species generator's 1e-11. All model/grid/time parameters are unchanged.
if nargin<1,error('Specify an output directory.');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
p=model.fig3_parameters;
% Run the missing three-species case first; each completed case is saved.
for m=[3,2]
 result=src.fig3_solve(p,m);
 save(fullfile(out_dir,sprintf('production_history_m%d.mat',m)),'result','-v7');
 d=result.diagnostics;
 vars={'t'};for j=1:m,vars{end+1}=sprintf('mass_fraction_%d',j);end
 writetable(array2table([d.time,d.mass_fraction],'VariableNames',vars), ...
  fullfile(out_dir,sprintf('production_mass_fractions_m%d.csv',m)));
 ix=find(abs(result.history.time-1)<1e-12,1);U=result.history.U(:,:,ix);
 vars={'x'};for j=1:m,vars{end+1}=sprintf('n%d',j);end;vars=[vars,{'n','p'}];
 writetable(array2table([result.x,U,sum(U,2),sum(U,2).^p.gamma],'VariableNames',vars), ...
  fullfile(out_dir,sprintf('production_profiles_t1_m%d.csv',m)));
 fprintf('Saved complete Figure3 production diagnostics and snapshots for m=%d.\n',m);
end
paths=utils.paths;template=fullfile(paths.data,'Figure3','archived_plot_data.mat');
if exist(template,'file'),utils.fig3_finalize_data(out_dir,template);end
end
