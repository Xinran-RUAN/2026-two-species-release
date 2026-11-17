function [par,plot_par]=two_dimensional_parameters(gamma,finalTime)
%TWO_DIMENSIONAL_PARAMETERS Original three-population model and grid.
par.name = 'Example 4: localized three-species expansion in 2D';
par.code_version = 'source-form-allocation-v6-2d-localized-clone-20260719';

% Domain and grid.
par.xmin = -10.0;
par.xmax =  10.0;
par.ymin = -10.0;
par.ymax =  10.0;
par.Nx = 500;
par.Ny = 500;

% Model parameters.
par.gamma = gamma;
par.pH = 1.0;
par.nH = par.pH^(1/par.gamma);
par.growth_rates = [8.0,4.0,1.5];
par.m = numel(par.growth_rates);

% Time discretization.
par.dt = 1.0e-2;
par.T = finalTime;

% Initial total density.
par.initial_radius = 2.0;
par.initial_density_height = 0.98;
par.initial_width = 0.08;

% Initial hidden composition.
par.fast_center = [2.05,0.00];
par.fast_sigma = 0.25;
par.fast_amplitude = 60.0;
par.slow_center = [-1.00,-0.90];
par.slow_sigma = 0.40;
par.slow_amplitude = 10.0;

% Only the snapshots used by the paper figure are stored.  Add further
% times here if more postprocessing data are needed.
par.snapshot_times = 0.0:0.1:par.T;

% Damped Newton parameters for the scalar density-pressure step.
par.newton_tol = 1.0e-10;
par.newton_maxit = 50;
par.line_search_maxit = 40;
par.armijo = 1.0e-4;
par.warn_residual_tol = 1.0e-8;
par.abort_on_newton_failure = true;

% Structure diagnostics and numerical-vacuum tolerances.
par.roundoff_tol = 1.0e-12;
par.energy_ineq_tol = 1.0e-7;
par.scalar_balance_tol = 1.0e-8;
par.remap_bound_tol = 1.0e-8;
par.remap_mass_tol = 1.0e-14;
par.remap_residual_tol = 1.0e-8;
par.reconstruction_tol = 1.0e-11;
par.occupied_tol = 1.0e-12;
par.abort_on_remap_failure = true;

% Progress and output.
par.progress_every = 10;
par.show_figures = true;
par.export_resolution = 300;

% Front and RGB-plot parameters.
par.front_level = 0.5;
par.front_angle_count = 361;
par.rgb_vacuum_level = 0.02;
par.rgb_opacity_power = 0.75;

% Plot parameters, collected here for easy adjustment.
plot_par.figure_position_rgb = [100,100,1350,430];
plot_par.figure_position_final = [120,120,1350,430];
plot_par.figure_position_diagnostics = [140,140,1100,430];
plot_par.font_size = 13;
plot_par.label_font_size = 15;
plot_par.title_font_size = 14;
plot_par.axes_line_width = 1.0;
plot_par.curve_line_width = 1.8;
plot_par.contour_line_width = 1.1;
plot_par.grid_alpha = 0.12;
plot_par.species_colors = [ ...
    0.70,0.05,0.32; ... % species 1: magenta-red
    0.95,0.58,0.08; ... % species 2: orange
    0.02,0.48,0.55];    % species 3: teal

end
