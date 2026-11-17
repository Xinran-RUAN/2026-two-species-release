function par=stiff_pressure_parameters
%STIFF_PRESSURE_PARAMETERS Fixed-grid seven-gamma comparison.
par.code_version = 'section4.4-current-scheme-flux-topology-20260727';

par.a = -5.0;
par.b =  5.0;
par.N = 6400;

par.dt = 1.0e-4;
par.T = 1.0;

par.pH = 1.0;
par.R0 = 1.0;
par.X0 = 0.25;

par.gammas = [10,20,40,80,160,320,640];

% Interior pressure residual:
% use cells for which p_infinity(x,T) >= interior_pressure_cutoff.
par.interior_pressure_cutoff = 0.2;

% Nonlinear scalar solve.
par.newton_tol = 1.0e-12;
par.newton_maxit = 40;
par.line_search_maxit = 25;
par.armijo = 1.0e-4;
par.warn_residual_tol = 1.0e-9;
par.abort_on_newton_failure = true;

% Roundoff and remapping diagnostics.
par.roundoff_tol = 1.0e-12;
par.remap_mass_tol = 1.0e-14;
par.remap_bound_tol = 1.0e-10;
par.remap_residual_tol = 1.0e-8;
par.reconstruction_tol = 1.0e-12;
par.abort_on_remap_failure = true;

% Output and progress.
par.resume = true;
par.progress_every = 500;
par.save_after_each_gamma = true;

end
