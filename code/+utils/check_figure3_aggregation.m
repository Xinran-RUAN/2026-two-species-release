function checks=check_figure3_aggregation
%CHECK_FIGURE3_AGGREGATION Small fixed-step consistency check, not a full run.
p=model.fig3_parameters;p.N=80;p.T=.001;p.snapshot_times=[0;p.T];p.progress_every=0;
a=src.fig3_solve(p,2);b=src.fig3_solve(p,3);
agg=[sum(b.final_U(:,1:2),2),b.final_U(:,3)];
checks.aggregation_max_abs_difference=max(abs(a.final_U(:)-agg(:)));
assert(checks.aggregation_max_abs_difference<5e-12);
q=p;q.a1=10;q.a2=2;q.store_solution_history=true;
init=model.structure_initial;
z=src.ms_tumor1d_solver(q,init);U=[z.n1,z.n2];
checks.legacy_two_species_max_abs_difference=max(abs(a.final_U(:)-U(:)));
assert(checks.legacy_two_species_max_abs_difference<5e-12);
checks.test_grid=p.N;checks.test_steps=round(p.T/p.dt);
checks.newton_tolerance=p.newton_tol;
disp(checks);
end
