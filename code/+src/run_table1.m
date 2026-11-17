function run_table1(out_dir)
%RUN_TABLE1 Repeat the paper fixed-grid timing benchmark from initial data.
% Three alternating timings follow one warm-up per scheme. Timesteps are the
% original accepted values, not selected by a new search. Wall-clock times
% depend on the current hardware, MATLAB version and system load.
if nargin<1,error('Specify an output directory for newly measured data.');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
cfg=model.paper_cases('Table1');
summary=struct([]);recorded_runs=struct([]);
for g=cfg.gammas
    p=model.accuracy_parameters('paper_smooth',g,cfg.N,cfg.reference_dt,cfg.T);
    ref=src.accuracy_solve(p,'dvd');assert(ref.success);
    p=model.accuracy_parameters('paper_smooth',g,cfg.N,cfg.warmup_dt,cfg.warmup_T);
    warm=src.accuracy_solve(p,'first');assert(warm.success);
    warm=src.accuracy_solve(p,'dvd');assert(warm.success);
    steps=cfg.first_order_steps(cfg.gammas==g);
    rows=struct([]);
    for rep=1:cfg.repetitions
        for methods={'first','dvd'}
            method=methods{1};p.T=cfg.T;
            if strcmp(method,'first'),p.dt=cfg.T/steps;else,p.dt=cfg.pc_dt;end
            o=src.accuracy_solve(p,method);assert(o.success);
            e=utils.s41_error(o.U,ref.U,o.h,g);r=o.stats;
            r.gamma=g;r.method=method;r.rep=rep;r.N=cfg.N;r.T=cfg.T;
            r.dt=p.dt;r.En=e.En;r.Esp=e.Esp;rows=append(rows,r);
            z=struct('gamma',g,'method',method,'repetition',rep,'out',o);
            recorded_runs=append(recorded_runs,z);
        end
    end
    writetable(struct2table(rows),fullfile(out_dir,sprintf('cost_g%d.csv',g)));
    L=rows(strcmp({rows.method},'first'));D=rows(strcmp({rows.method},'dvd'));
    r=struct('gamma',g,'N',cfg.N,'T',cfg.T,'dt_first',cfg.T/steps,'dt_dvd',cfg.pc_dt, ...
        'Esp_first',L(1).Esp,'Esp_dvd',D(1).Esp, ...
        'seconds_first',median([L.seconds]),'seconds_dvd',median([D.seconds]), ...
        'speedup',median([L.seconds])/median([D.seconds]),'repetitions',cfg.repetitions);
    summary=append(summary,r);
    save(fullfile(out_dir,'table1_recomputed_runs.mat'),'recorded_runs','summary','-v7');
end
writetable(struct2table(summary),fullfile(out_dir,'cost.csv'));
fprintf('Finished new Table 1 timing measurements. Historical timings were not overwritten.\n');
end
function v=append(v,r)
if isempty(v),v=r;else,v(end+1)=r;end
end
