function p=fig3_parameters
%FIG3_PARAMETERS Parameters from the final submitted two/three-species test.
p.a=-10;p.b=10;p.N=2000;p.gamma=20;p.pH=1;p.nH=1;p.T=2;p.dt=1e-4;
p.newton_tol=1e-12;p.newton_maxit=50;p.line_search_maxit=30;p.armijo=1e-4;
p.warn_residual_tol=1e-8;p.abort_on_newton_failure=true;
p.roundoff_tol=1e-12;p.remap_mass_tol=1e-13;p.remap_negative_tol=1e-10;p.occupied_tol=1e-12;
p.snapshot_times=(0:.1:2)';p.progress_every=500;
end
