function run_figure1(out_dir)
%RUN_FIGURE1 Recompute the paper accuracy experiment from its initial data.
% This is an explicitly requested long computation, not a redraw operation.
% All timesteps, grids, initial cell averages and solver settings agree with
% the original run_convergence_experiments/run_time_error_tests workflows.
if nargin<1,error('Specify an output directory for newly computed data.');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
cfg=model.paper_cases('Figure1');
cases={'smooth','vacuum'};temporal=struct([]);spatial=struct([]);
refs=cell(2,1);test_states=cell(2,2,7);spatial_states=cell(2,5);
dts=cfg.time_steps;Ns=cfg.spatial_grids;
for c=1:2
    kind='paper_smooth';if c==2,kind='paper_structure';end
    p=model.accuracy_parameters(kind,cfg.gamma,cfg.temporal_cells,cfg.reference_dt,cfg.final_time);
    p.initial_species=model.fallback_initial(p);
    refs{c}=src.accuracy_solve(p,'dvd');assert(refs{c}.success);
    for im=1:2
        methods={'first','dvd'};method=methods{im};previous=NaN;
        for j=1:numel(dts)
            p.dt=dts(j);o=src.accuracy_solve(p,method);assert(o.success);
            test_states{c,im,j}=o;e=utils.s41_error(o.U,refs{c}.U,o.h,cfg.gamma);
            order=NaN;if j>1,order=log2(previous/e.Esp);end
            r=struct('case_name',cases{c},'method',method,'dt',p.dt, ...
                'Esp',e.Esp,'observed_order',order, ...
                'fallback_count',o.stats.fallback_steps);
            temporal=append(temporal,r);previous=e.Esp;
        end
    end
    for j=1:numel(Ns)
        p=model.accuracy_parameters(kind,cfg.gamma,Ns(j),cfg.spatial_dt,cfg.final_time);
        p.initial_species=model.fallback_initial(p);
        spatial_states{c,j}=src.accuracy_solve(p,'first');assert(spatial_states{c,j}.success);
    end
    previous=NaN;
    for j=1:numel(Ns)-1
        o=spatial_states{c,j};fine=spatial_states{c,j+1}.U;
        V=(fine(1:2:end,:)+fine(2:2:end,:))/2;e=utils.s41_error(o.U,V,o.h,cfg.gamma);
        order=NaN;if j>1,order=log2(previous/e.Esp);end
        r=struct('case_name',cases{c},'N',Ns(j),'Nfine',Ns(j+1), ...
            'h',o.h,'dt',cfg.spatial_dt,'D_sp',e.L1_n1+e.L1_n2, ...
            'E_sp_space',e.Esp,'observed_order',order);
        spatial=append(spatial,r);previous=e.Esp;
    end
    save(fullfile(out_dir,'figure1_recomputed_states.mat'), ...
        'refs','test_states','spatial_states','temporal','spatial','dts','Ns','-v7');
end
writetable(struct2table(temporal),fullfile(out_dir,'temporal_errors.csv'));
writetable(struct2table(spatial),fullfile(out_dir,'spatial_errors.csv'));
fprintf('Finished Figure 1 fixed-step computations. Published data were not overwritten.\n');
end
function v=append(v,r)
if isempty(v),v=r;else,v(end+1)=r;end
end
