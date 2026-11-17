function par=switching_parameters
%SWITCHING_PARAMETERS Pressure-dependent phenotype-switching test.
par.name = 'Example 2: pressure-induced phenotype switching';
par.code_version = 'source-form-allocation-v6-transition-20260727';

par.a = -10.0;
par.b =  10.0;
par.N = 2000;

par.gamma = 20;
par.pH = 1.0;
par.nH = par.pH^(1/par.gamma);

par.T = 2.0;
par.dt = 1.0e-4;

% Source-wise net growth rates:
%   H_1(p)=H_2(p)=a_H(p_H-p).
par.aH = 6.0;

% Switching rates:
%   type 1 -> type 2: alpha*p,
%   type 2 -> type 1: beta*(p_H-p).
par.alpha = 1.0;
par.beta = 0.25;

% At p=p_H/2, the frozen switching equilibrium is
% beta/(alpha+beta).  This is only a reference value; the equilibrium is
% pressure dependent.
par.c_star = par.beta/(par.alpha+par.beta);

% Use the same total-density profile as the structure/non-closure test.
par.ell_n = 0.05;

% Dense snapshots are retained for the space-time composition plot.
par.snapshot_times = 0:0.01:par.T;

% Damped Newton parameters for the density-pressure step.
par.newton_tol = 1.0e-12;
par.newton_maxit = 50;
par.line_search_maxit = 30;
par.armijo = 1.0e-4;
par.warn_residual_tol = 1.0e-8;
par.abort_on_newton_failure = true;

% Structure diagnostics.
par.roundoff_tol = 1.0e-12;
par.energy_ineq_tol = 1.0e-8;
par.scalar_balance_tol = 1.0e-8;
par.remap_bound_tol = 1.0e-8;
par.remap_mass_tol = 1.0e-14;
par.remap_residual_tol = 1.0e-8;
par.reconstruction_tol = 1.0e-12;
par.abort_on_remap_failure = true;
par.occupied_tol = 0.0;
par.progress_every = 100;

% Diagnostics for the high-pressure core and low-pressure front.
par.core_pressure_threshold = 0.8;
par.front_pressure_low = 0.05;
par.front_pressure_high = 0.3;

% Plot options for the single paper heat map.
par.show_figures = true;
par.plot_species = 2;  % paper figure displays phi_2
par.vacuum_tol = 1.0e-4;

par.font_size = 14;
par.label_font_size = 16;
par.title_font_size = 16;
par.legend_font_size = 11;
par.axis_line_width = 1.0;

par.figure_position = [100,100,760,500];
par.heatmap_axes_position = [0.105,0.145,0.700,0.775];
par.colorbar_position = [0.845,0.145,0.026,0.775];

par.plot_xlim = [-8,8];

% phi_2 is a fraction, so [0,1] is the natural physical range.
% Keeping the full range also prevents depletion below 0.5 from being
% hidden by color saturation.
par.phi_color_limits = [0,1];
par.phi_colorbar_ticks = 0:0.2:1;

% Pressure-contour styles.
par.core_contour_color = [0,0,0];
par.core_contour_style = '-';
par.core_contour_width = 1.5;

par.front_high_contour_color = [0.20,0.20,0.20];
par.front_high_contour_style = '--';
par.front_high_contour_width = 1.3;

par.front_low_contour_color = [0.20,0.20,0.20];
par.front_low_contour_style = ':';
par.front_low_contour_width = 1.3;

par.export_resolution = 300;

end
