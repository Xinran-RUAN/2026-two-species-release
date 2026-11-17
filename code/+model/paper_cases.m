function cfg=paper_cases(case_id)
%PAPER_CASES Production cases and unchanged experiment settings.
% MODEL.PAPER_CASES lists the supported initial-value computations.
ids={'Figure1','Figure2','Figure3','Table1','Figure4','Figure5','Figure6Gamma3'};
if nargin<1,cfg=ids;return;end
index=find(strcmpi(char(case_id),ids),1);
assert(~isempty(index),'two_species:UnknownCase','Unknown production case: %s',char(case_id));
cfg.id=ids{index};cfg.output_folder=cfg.id;
switch cfg.id
 case 'Figure1'
  cfg.gamma=20;cfg.temporal_cells=1000;cfg.final_time=.2;cfg.reference_dt=5e-6;
  cfg.time_steps=[.008,.004,.002,.001,.0005,.00025,.000125];
  cfg.spatial_grids=[100,200,400,800,1600];cfg.spatial_dt=1e-5;
 case 'Figure2'
  cfg.parameters=model.structure_parameters;
 case 'Figure3'
  cfg.parameters=model.fig3_parameters;
 case 'Table1'
  cfg.gammas=[3,20];cfg.N=1000;cfg.T=.2;cfg.reference_dt=1/25600;
  cfg.pc_dt=.01;cfg.first_order_steps=[5123,569];cfg.repetitions=3;
  cfg.warmup_dt=.001;cfg.warmup_T=.01;
 case 'Figure4'
  cfg.parameters=model.switching_parameters;
 case 'Figure5'
  cfg.parameters=model.stiff_pressure_parameters;
  cfg.additional_outputs={'Table2','Table3'};
 case 'Figure6Gamma3'
  cfg.gamma=3;cfg.final_time=2.5;cfg.output_folder='Figure6';
  [cfg.parameters,cfg.plot_parameters]=model.two_dimensional_parameters(cfg.gamma,cfg.final_time);
end
end
