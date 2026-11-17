function par=structure_parameters
%STRUCTURE_PARAMETERS Archived first-order structure-test parameters.
par.a=-10;par.b=10;par.N=2000;par.gamma=20;par.pH=1;par.nH=1;
par.a1=10;par.a2=2;par.T=2;par.dt=1e-4;
par.newton_tol=1e-11;par.newton_maxit=50;par.warn_residual_tol=1e-8;
par.progress_every=500;par.store_solution_history=true;
snapshot_dt=.1;par.snapshot_times=unique([(0:snapshot_dt:par.T)';par.T]);
par.occupied_tol=1e-12;par.remap_mass_tol=1e-13;par.remap_negative_tol=1e-10;
end
