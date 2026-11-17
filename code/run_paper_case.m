function case_dir=run_paper_case(case_id,runs_root)
%RUN_PAPER_CASE Recompute one paper case into a separate numerical data tree.
% RUN_PAPER_CASE('Figure3',RUNS) writes RUNS/Figure3. The resulting data tree
% can be used by REDRAW_PAPER('Figure3',OUT,RUNS).
% MODEL.PAPER_CASES lists the supported computations.
paths=startup_paper;
if nargin<1,error('two_species:MissingCase','Select a case from model.paper_cases.');end
if nargin<2 || isempty(runs_root),runs_root=fullfile(tempdir,'two_species_paper_runs');end
cfg=model.paper_cases(case_id);runs_root=char(runs_root);
a=char(java.io.File(paths.data).getCanonicalPath());b=char(java.io.File(runs_root).getCanonicalPath());
assert(~strcmp(a,b) && ~startsWith(b,[a filesep]), ...
 'two_species:ReadOnlyData','Select a run directory outside the published data tree.');
case_dir=fullfile(runs_root,cfg.output_folder);
if ~isfolder(case_dir),mkdir(case_dir);end
switch cfg.id
 case 'Figure1'
  src.run_figure1(case_dir);
  utils.finalize_accuracy_case(case_dir,fullfile(paths.data,'Figure1','plot_data.mat'));
 case 'Figure2'
  src.run_figure2(case_dir);
  utils.finalize_structure_case(case_dir,fullfile(paths.data,'Figure2','plot_data.mat'));
 case 'Figure3'
  src.run_figure3(case_dir);
 case 'Table1'
  src.run_table1(case_dir);
 case 'Figure4'
  src.run_example2_transition_1d(case_dir);
  copyfile(fullfile(case_dir,'example2_data_conservative_remap.mat'),fullfile(case_dir,'switching_data.mat'));
 case 'Figure5'
  src.run_section4_4_stiff_pressure(case_dir);
  source=fullfile(case_dir,'section4_4_stiff_pressure_data.mat');
  copyfile(source,fullfile(case_dir,'stiff_pressure_data.mat'));
  for id={'Table2','Table3'}
   dest=fullfile(runs_root,id{1});if ~isfolder(dest),mkdir(dest);end
   copyfile(source,fullfile(dest,'stiff_pressure_data.mat'));
  end
 case 'Figure6Gamma3'
  case_dir=fullfile(case_dir,'gamma3');
  src.run_example4_2d_multispecies_v5(cfg.gamma,cfg.final_time,case_dir);
  copyfile(fullfile(case_dir,'example4_2d_multispecies_data.mat'), ...
   fullfile(runs_root,'Figure6','gamma3_data.mat'));
end
fprintf('Computed %s in %s\n',cfg.id,case_dir);
end
