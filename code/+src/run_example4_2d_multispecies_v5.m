function run_example4_2d_multispecies_v5(gamma,finalTime,outdir)
%RUN_EXAMPLE4_2D_MULTISPECIES_V5
% Example 4: localized three-species expansion in two dimensions.
%
% Numerical scheme used here:
%   1. scalar Onsager density-pressure step with time-lagged mobility;
%   2. local normalized exact frozen-pressure reaction allocation;
%   3. implicit conservative upwind species remapping driven by exactly
%      the same physical total-mass flux as the density step.
%
% The fully discrete state is the species-density array
%   n_species(:,:,s),  s=1,...,m.
% The total density and pressure are reconstructed from
%   n = sum_s n_s,  p=n^gamma.
% Composition fractions are derived only in numerically occupied cells and
% are never transported as independent all-grid unknowns.  This is the same
% conservative-remapping scheme used in the corrected Example 1 and in the
% current paper.
%
% Default test:
%   m=3,
%   H_s(p)=a_s(p_H-p),  a=[8,4,1.5],
% with no phenotype conversion.  A localized fast-growing clone is placed
% near the initial right front, a medium-growth clone forms the background,
% and a slower clone is placed in the lower-left interior.
%
% Use run_paper_case('Figure6Gamma3',runs) for the paper case.
% Output is written to the directory supplied by the caller.

clc; close all;
% run_paper_case selects the paper final time.
if nargin < 1 || isempty(gamma), gamma=3; end
if nargin < 2 || isempty(finalTime), finalTime=3; end
if nargin < 3 || isempty(outdir), outdir=fullfile(tempdir,'two_species_rerun',sprintf('2d_gamma%g',gamma)); end

%% ========================================================================
%  Parameters
%  ========================================================================
[par,plot_par]=model.two_dimensional_parameters(gamma,finalTime);
fprintf('Code version: %s\n',par.code_version);

%% ========================================================================
%  Output folder
%  ========================================================================
% outdir is supplied by the caller.
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% ========================================================================
%  Grid and time
%  ========================================================================
Nx = par.Nx;
Ny = par.Ny;
hx = (par.xmax-par.xmin)/Nx;
hy = (par.ymax-par.ymin)/Ny;
if abs(hx-hy) > 1.0e-14*max(1,max(abs([hx,hy])))
    error('This implementation assumes hx=hy.');
end
h = hx;

x = par.xmin+((1:Nx)-0.5)*h;
y = par.ymin+((1:Ny)-0.5)*h;
[X,Y] = meshgrid(x,y);

nt = round(par.T/par.dt);
par.dt = par.T/nt;
time_grid = (0:nt)'*par.dt;
newton_time_grid = time_grid(2:end);

snapshot_steps = round(par.snapshot_times/par.dt);
snapshot_steps = max(0,min(nt,snapshot_steps));
[snapshot_steps,ia] = unique(snapshot_steps,'stable');
snapshot_times = snapshot_steps(:)'*par.dt;
par.snapshot_times = par.snapshot_times(ia);
nsnapshot = numel(snapshot_steps);

assert(~isempty(snapshot_steps) && snapshot_steps(1)==0, ...
    'The snapshot list must include t=0.');
assert(snapshot_steps(end)==nt, ...
    'The snapshot list must include the final time.');

%% ========================================================================
%  Initial total density and species densities
%  ========================================================================
[n0,phi0]=model.two_dimensional_initial(X,Y,par);
phi0 = normalize_simplex_array_2d(phi0,par.roundoff_tol);

n_species = zeros(Ny,Nx,par.m);
for s = 1:par.m
    n_species(:,:,s) = n0.*phi0(:,:,s);
end

n = sum(n_species,3);
p = n.^par.gamma;

initial_reconstruction_error = max(abs(n(:)-n0(:)));
if initial_reconstruction_error > par.reconstruction_tol
    error('Initial species reconstruction error is %.3e.', ...
        initial_reconstruction_error);
end

%% ========================================================================
%  Snapshot storage
%  ========================================================================
snaps = struct();
snaps.t = zeros(1,nsnapshot);
snaps.n = cell(1,nsnapshot);
snaps.p = cell(1,nsnapshot);
snaps.n_species = cell(1,nsnapshot);

snapshot_id = 1;
if snapshot_steps(snapshot_id)==0
    snaps = store_snapshot_2d(snaps,snapshot_id,0,n,p,n_species);
    snapshot_id = snapshot_id+1;
end

%% ========================================================================
%  Diagnostics
%  ========================================================================
energy_hist = zeros(nt+1,1);
mass_hist = zeros(nt+1,1);
species_mass_hist = zeros(nt+1,par.m);
min_species_hist = zeros(nt+1,1);
min_n_hist = zeros(nt+1,1);
max_n_hist = zeros(nt+1,1);
max_p_hist = zeros(nt+1,1);
reconstruction_err_hist = zeros(nt+1,1);
simplex_err_hist = zeros(nt+1,1);

newton_it_hist = zeros(nt,1);
newton_res_hist = zeros(nt,1);
transport_diss_hist = zeros(nt,1);
reaction_diss_hist = zeros(nt,1);
energy_ineq_defect_hist = zeros(nt,1);
source_balance_defect_hist = zeros(nt,1);
reaction_allocation_err_hist = zeros(nt,1);
remap_simplex_err_hist = zeros(nt,1);
remap_reconstruction_err_hist = zeros(nt,1);
remap_residual_hist = zeros(nt,1);
remap_mass_compatibility_hist = zeros(nt,1);
active_cell_count_hist = zeros(nt,1);
zero_mass_fallback_count_hist = zeros(nt,1);
max_raw_negative_rhs_hist = zeros(nt,1);

remap_overshoot_count = 0;
remap_overshoot_max = 0.0;

[state_diag,species_masses] = state_diagnostics_2d( ...
    n_species,n,p,h,par);
energy_hist(1) = modified_energy_2d(n,h,par.gamma,par.pH);
mass_hist(1) = h^2*sum(n(:));
species_mass_hist(1,:) = species_masses;
min_species_hist(1) = state_diag.min_species;
min_n_hist(1) = state_diag.min_n;
max_n_hist(1) = state_diag.max_n;
max_p_hist(1) = state_diag.max_p;
reconstruction_err_hist(1) = state_diag.reconstruction_error;
simplex_err_hist(1) = state_diag.simplex_error;

fprintf(['Example 4: 2D %d-species conservative-remapping run, ', ...
    'Nx=%d, Ny=%d, steps=%d.\n'],par.m,Nx,Ny,nt);

%% ========================================================================
%  Time stepping
%  ========================================================================
for k = 1:nt
    n_species_old = n_species;
    nold = sum(n_species_old,3);
    Eold = energy_hist(k);

    old_reconstruction_array = abs(sum(n_species_old,3)-nold);
    old_reconstruction_error = max(old_reconstruction_array(:));
    if old_reconstruction_error > par.reconstruction_tol
        error('Reconstruction error before step %d is %.3e.', ...
            k,old_reconstruction_error);
    end

    % Frozen scalar reaction mobility:
    %   M^k=sum_s a_s n_s^k.
    Mmob = reaction_mobility_species_2d(n_species_old, ...
        par.growth_rates);

    % Step 1: scalar Onsager density-pressure update.
    [nnew,pnew,newton_info] = solve_density_step_2d( ...
        nold,Mmob,par,h);
    newton_it_hist(k) = newton_info.iter;
    newton_res_hist(k) = newton_info.residual;

    Enew = modified_energy_2d(nnew,h,par.gamma,par.pH);
    [Dtrans,Dreact,Edefect] = density_energy_diagnostics_2d( ...
        nold,pnew,Mmob,Eold,Enew,par,h);
    transport_diss_hist(k) = Dtrans;
    reaction_diss_hist(k) = Dreact;
    energy_ineq_defect_hist(k) = Edefect;

    if Edefect > par.energy_ineq_tol
        warning('Energy inequality defect at step %d is %.3e.', ...
            k,Edefect);
    end

    % Physical total-mass flux J=-n^k grad_h p^{k+1}.
    [Jx,Jy] = physical_mass_flux_2d(nold,pnew,h);

    % In an exact scalar solve, the source and conservative forms agree.
    S = Mmob.*(par.pH-pnew);
    source_scale = max(1,max(abs(S(:))));
    source_tol = par.roundoff_tol*source_scale;
    if any(S(:)<-source_tol)
        error('Negative scalar source at step %d: min(S)=%.3e.', ...
            k,min(S(:)));
    end
    S(S<0) = 0;

    nhat_source = nold+par.dt*S;
    divJ = divergence_flux_2d(Jx,Jy,h);
    nhat_balance = nnew+par.dt*divJ;
    source_balance_defect = max(abs(nhat_balance(:)-nhat_source(:)));
    source_balance_defect_hist(k) = source_balance_defect;

    adaptive_balance_tol = max( ...
        par.scalar_balance_tol, ...
        5.0*newton_info.residual+100*eps);
    if source_balance_defect > adaptive_balance_tol
        warning(['Scalar balance defect at step %d is %.3e ', ...
            '(Newton residual %.3e).'], ...
            k,source_balance_defect,newton_info.residual);
    end

    % Reaction allocation uses the analytically nonnegative source form;
    % the conservative form above is retained only as a scalar diagnostic.
    nhat = nhat_source;

    nhat_scale = max(1.0,max(abs(nhat(:))));
    nhat_tol = max( ...
        par.roundoff_tol*nhat_scale, ...
        100*eps(nhat_scale)*nhat_scale);
    small_negative = nhat < 0 & nhat >= -nhat_tol;
    nhat(small_negative) = 0;

    if any(nhat(:) < -nhat_tol)
        error(['Physical source-updated mass is negative at ', ...
            'step %d: min=%.3e.'], ...
            k,min(nhat(:)));
    end

    % Step 2: local normalized exact frozen-pressure reaction allocation.
    [nhat_species,reaction_info] = ...
        local_reaction_allocation_diagonal_2d( ...
        n_species_old,nhat,pnew,par);
    reaction_allocation_err_hist(k) = ...
        reaction_info.allocation_error;

    % Step 3: implicit conservative upwind species remapping using J.
    [n_species,theta,remap_info] = ...
        conservative_species_remap_2d( ...
        nnew,nhat_species,Jx,Jy,pnew,par.dt,h,par);

    n = nnew;
    p = pnew;

    remap_simplex_err_hist(k) = remap_info.simplex_error;
    remap_reconstruction_err_hist(k) = ...
        remap_info.reconstruction_error;
    remap_residual_hist(k) = remap_info.species_residual;
    remap_mass_compatibility_hist(k) = ...
        remap_info.mass_compatibility_defect;
    active_cell_count_hist(k) = remap_info.active_cell_count;
    zero_mass_fallback_count_hist(k) = ...
        remap_info.zero_mass_fallback_count;
    max_raw_negative_rhs_hist(k) = remap_info.max_raw_negative_rhs;

    if remap_info.overshoot > par.remap_bound_tol
        remap_overshoot_count = remap_overshoot_count+1;
    end
    remap_overshoot_max = max(remap_overshoot_max, ...
        remap_info.overshoot);

    [state_diag,species_masses] = state_diagnostics_2d( ...
        n_species,n,p,h,par);
    energy_hist(k+1) = Enew;
    mass_hist(k+1) = h^2*sum(n(:));
    species_mass_hist(k+1,:) = species_masses;
    min_species_hist(k+1) = state_diag.min_species;
    min_n_hist(k+1) = state_diag.min_n;
    max_n_hist(k+1) = state_diag.max_n;
    max_p_hist(k+1) = state_diag.max_p;
    reconstruction_err_hist(k+1) = ...
        state_diag.reconstruction_error;
    simplex_err_hist(k+1) = state_diag.simplex_error;

    if snapshot_id<=nsnapshot && k==snapshot_steps(snapshot_id)
        snaps = store_snapshot_2d(snaps,snapshot_id,k*par.dt, ...
            n,p,n_species);
        snapshot_id = snapshot_id+1;
    end

    if mod(k,par.progress_every)==0 || k==1 || k==nt
        fprintf([' step %4d/%4d, t=%.3f, Newton=%2d, ', ...
            'res=%.3e, max p=%.8f, Edef=%.3e, ', ...
            'remap=%.3e, recon=%.3e\n'], ...
            k,nt,k*par.dt,newton_info.iter,newton_info.residual, ...
            max_p_hist(k+1),Edefect,remap_info.species_residual, ...
            remap_info.reconstruction_error);
    end
end

if snapshot_id~=nsnapshot+1
    error('Not every requested snapshot was saved.');
end

%% ========================================================================
%  Final fields and morphology diagnostics
%  ========================================================================
phi = species_fractions_2d(n_species,n,par.occupied_tol);

[front_theta,front_radius_final] = compute_front_radius_2d( ...
    x,y,n,par.front_level,par.front_angle_count);

valid_front = isfinite(front_radius_final);
if any(valid_front)
    max_front_radius = max(front_radius_final(valid_front));
    min_front_radius = min(front_radius_final(valid_front));
    front_radius_range = max_front_radius-min_front_radius;
    protrusion_amplitude = max_front_radius- ...
        mean(front_radius_final(valid_front));
else
    max_front_radius = NaN;
    min_front_radius = NaN;
    front_radius_range = NaN;
    protrusion_amplitude = NaN;
end

%% ========================================================================
%  Summary and pass/fail criteria
%  ========================================================================
summary.min_species = min(min_species_hist);
summary.min_n = min(min_n_hist);
summary.max_n = max(max_n_hist);
summary.max_p = max(max_p_hist);
summary.max_pressure_overshoot = max(0,summary.max_p-par.pH);
summary.max_reconstruction_error = max(reconstruction_err_hist);
summary.max_simplex_error = max(simplex_err_hist);
summary.max_positive_energy_jump = max([0;diff(energy_hist)]);
summary.max_energy_inequality_defect = max(energy_ineq_defect_hist);
summary.max_scalar_balance_defect = max(source_balance_defect_hist);
summary.max_reaction_allocation_error = ...
    max(reaction_allocation_err_hist);
summary.max_remap_simplex_error = max(remap_simplex_err_hist);
summary.max_remap_reconstruction_error = ...
    max(remap_reconstruction_err_hist);
summary.max_remap_residual = max(remap_residual_hist);
summary.max_remap_mass_compatibility = ...
    max(remap_mass_compatibility_hist);
summary.max_newton_residual = max(newton_res_hist);
summary.mean_newton_iterations = mean(newton_it_hist);
summary.max_newton_iterations = max(newton_it_hist);
summary.remap_overshoot_count = remap_overshoot_count;
summary.remap_overshoot_max = remap_overshoot_max;
summary.initial_energy = energy_hist(1);
summary.final_energy = energy_hist(end);
summary.final_energy_drop = energy_hist(end)-energy_hist(1);
summary.final_total_mass = mass_hist(end);
summary.max_front_radius = max_front_radius;
summary.min_front_radius = min_front_radius;
summary.front_radius_range = front_radius_range;
summary.protrusion_amplitude = protrusion_amplitude;

criteria.species_nonnegative = ...
    summary.min_species>=-par.roundoff_tol;
criteria.density_bounds = ...
    summary.min_n>=-par.roundoff_tol && ...
    summary.max_n<=par.nH+1.0e-10;
criteria.pressure_bound = summary.max_p<=par.pH+1.0e-10;
criteria.reconstruction = ...
    summary.max_reconstruction_error<=par.reconstruction_tol;
criteria.simplex = summary.max_simplex_error<=par.remap_bound_tol;
criteria.energy_inequality = ...
    summary.max_energy_inequality_defect<=par.energy_ineq_tol && ...
    summary.max_positive_energy_jump<=par.energy_ineq_tol;
criteria.scalar_balance = ...
    summary.max_scalar_balance_defect<=par.scalar_balance_tol;
criteria.reaction_allocation = ...
    summary.max_reaction_allocation_error<=par.reconstruction_tol;
criteria.remap_residual = ...
    summary.max_remap_residual<=par.remap_residual_tol;
criteria.newton_residual = ...
    summary.max_newton_residual<=par.warn_residual_tol;
criteria.all_pass = ...
    criteria.species_nonnegative && ...
    criteria.density_bounds && ...
    criteria.pressure_bound && ...
    criteria.reconstruction && ...
    criteria.simplex && ...
    criteria.energy_inequality && ...
    criteria.scalar_balance && ...
    criteria.reaction_allocation && ...
    criteria.remap_residual && ...
    criteria.newton_residual;

%% ========================================================================
%  Save data
%  ========================================================================
diagnostics = struct();
diagnostics.time = time_grid;
diagnostics.energy = energy_hist;
diagnostics.mass = mass_hist;
diagnostics.species_mass = species_mass_hist;
diagnostics.min_species = min_species_hist;
diagnostics.min_n = min_n_hist;
diagnostics.max_n = max_n_hist;
diagnostics.max_p = max_p_hist;
diagnostics.reconstruction_error = reconstruction_err_hist;
diagnostics.simplex_error = simplex_err_hist;
diagnostics.newton_time = newton_time_grid;
diagnostics.newton_iterations = newton_it_hist;
diagnostics.newton_residual = newton_res_hist;
diagnostics.transport_dissipation = transport_diss_hist;
diagnostics.reaction_dissipation = reaction_diss_hist;
diagnostics.energy_inequality_defect = energy_ineq_defect_hist;
diagnostics.scalar_balance_defect = source_balance_defect_hist;
diagnostics.reaction_allocation_error = reaction_allocation_err_hist;
diagnostics.remap_simplex_error = remap_simplex_err_hist;
diagnostics.remap_reconstruction_error = ...
    remap_reconstruction_err_hist;
diagnostics.remap_residual = remap_residual_hist;
diagnostics.remap_mass_compatibility = ...
    remap_mass_compatibility_hist;
diagnostics.active_cell_count = active_cell_count_hist;
diagnostics.zero_mass_fallback_count = zero_mass_fallback_count_hist;
diagnostics.max_raw_negative_rhs = max_raw_negative_rhs_hist;
diagnostics.front_theta = front_theta;
diagnostics.front_radius_final = front_radius_final;
diagnostics.front_level = par.front_level;

params = par;
data_file = fullfile(outdir,'example4_2d_multispecies_data.mat');
save(data_file, ...
    'x','y','X','Y','n0','phi0','n','p','n_species','phi', ...
    'snaps','diagnostics','summary','criteria','params', ...
    'front_theta','front_radius_final','-v7.3');

summary_table = struct2table(summary);
criteria_table = struct2table(criteria);
writetable(summary_table,fullfile(outdir,'example4_summary.csv'));
writetable(criteria_table,fullfile(outdir,'example4_pass_fail.csv'));
utils.write_time_diagnostics_csv_2d( ...
    fullfile(outdir,'example4_time_diagnostics.csv'), ...
    time_grid,energy_hist,max_p_hist,reconstruction_err_hist, ...
    simplex_err_hist,newton_it_hist,newton_res_hist, ...
    energy_ineq_defect_hist,source_balance_defect_hist, ...
    remap_residual_hist);

%% ========================================================================
%  Plots
%  ========================================================================
%% ========================================================================
%  Command-line summary
%  ========================================================================
fprintf('\n================ EXAMPLE 4 SUMMARY ================\n');
fprintf('min species density              = %.12e\n',summary.min_species);
fprintf('min/max total density            = %.12e / %.12e\n', ...
    summary.min_n,summary.max_n);
fprintf('max pressure                     = %.12e\n',summary.max_p);
fprintf('max pressure overshoot           = %.12e\n', ...
    summary.max_pressure_overshoot);
fprintf('max reconstruction error         = %.12e\n', ...
    summary.max_reconstruction_error);
fprintf('max occupied simplex error       = %.12e\n', ...
    summary.max_simplex_error);
fprintf('max positive energy jump         = %.12e\n', ...
    summary.max_positive_energy_jump);
fprintf('max energy inequality defect     = %.12e\n', ...
    summary.max_energy_inequality_defect);
fprintf('max scalar balance defect        = %.12e\n', ...
    summary.max_scalar_balance_defect);
fprintf('max reaction allocation error    = %.12e\n', ...
    summary.max_reaction_allocation_error);
fprintf('max remapping residual           = %.12e\n', ...
    summary.max_remap_residual);
fprintf('max remap compatibility defect   = %.12e\n', ...
    summary.max_remap_mass_compatibility);
fprintf('max Newton residual              = %.12e\n', ...
    summary.max_newton_residual);
fprintf('Newton iterations mean/max       = %.2f / %d\n', ...
    summary.mean_newton_iterations,summary.max_newton_iterations);
fprintf('front radius min/max/range       = %.6f / %.6f / %.6f\n', ...
    min_front_radius,max_front_radius,front_radius_range);
fprintf('protrusion amplitude             = %.6f\n', ...
    protrusion_amplitude);
fprintf('OVERALL PASS                     = %d\n',criteria.all_pass);
fprintf('Data file: %s\n',data_file);
fprintf('Output folder: %s\n',outdir);
fprintf('===================================================\n');

end

%% ========================================================================
function [nnew,pnew,info] = solve_density_step_2d(nold,Mmob,par,h)
%SOLVE_DENSITY_STEP_2D Projected damped Newton solve of
%
%   n^{k+1}-n^k
%   = dt div_h(n^k_edge grad_h p^{k+1})
%     -dt M^k(p^{k+1}-p_H),
%   p^{k+1}=(n^{k+1})^gamma.

[Ny,Nx] = size(nold);
n = min(max(nold,0),par.nH);

[Ax,Ay] = density_edge_mobilities_2d(nold);

info.converged = false;
info.iter = 0;
info.residual = Inf;
info.total_line_search_steps = 0;
info.min_accepted_lambda = 1.0;

stop_tol = par.newton_tol*(1+norm(nold(:),inf));

for it = 1:par.newton_maxit
    [R,Jmat] = density_residual_jacobian_2d( ...
        n,nold,Ax,Ay,Mmob,par,h);
    res_inf = norm(R(:),inf);

    if res_inf<=stop_tol
        info.converged = true;
        info.iter = it-1;
        info.residual = res_inf;
        break;
    end

    dn_vec = -(Jmat\R(:));
    if any(~isfinite(dn_vec))
        error('Density Newton produced a nonfinite search direction.');
    end
    dn = reshape(dn_vec,Ny,Nx);

    lambda = 1.0;
    accepted = false;
    res2 = norm(R(:),2);
    best_n = n;
    best_res2 = res2;
    best_ls = 0;

    for ls = 1:par.line_search_maxit
        ntrial = n+lambda*dn;
        ntrial = min(max(ntrial,0),par.nH);

        Rtrial = density_residual_2d( ...
            ntrial,nold,Ax,Ay,Mmob,par,h);
        trial_res2 = norm(Rtrial(:),2);

        if isfinite(trial_res2) && trial_res2<best_res2
            best_res2 = trial_res2;
            best_n = ntrial;
            best_ls = ls;
        end

        if isfinite(trial_res2) && ...
                trial_res2<=(1-par.armijo*lambda)*res2
            accepted = true;
            break;
        end
        lambda = 0.5*lambda;
    end

    info.total_line_search_steps = ...
        info.total_line_search_steps+ls;

    if accepted
        n = ntrial;
        info.iter = it;
        info.residual = norm(Rtrial(:),inf);
        info.min_accepted_lambda = ...
            min(info.min_accepted_lambda,lambda);
    elseif best_res2<res2
        n = best_n;
        info.iter = it;
        Rbest = density_residual_2d( ...
            n,nold,Ax,Ay,Mmob,par,h);
        info.residual = norm(Rbest(:),inf);
        if best_ls>0
            info.min_accepted_lambda = min( ...
                info.min_accepted_lambda,2^(-(best_ls-1)));
        end
    else
        info.iter = it;
        info.residual = res_inf;
        warning(['Density Newton line search failed at iteration %d; ', ...
            'residual=%.3e.'],it,res_inf);
        break;
    end
end

nnew = roundoff_project_interval_2d( ...
    n,0,par.nH,par.roundoff_tol,'density n');
pnew = nnew.^par.gamma;

Rfinal = density_residual_2d( ...
    nnew,nold,Ax,Ay,Mmob,par,h);
info.residual = norm(Rfinal(:),inf);
if info.residual<=stop_tol
    info.converged = true;
end

if ~info.converged && info.residual>par.warn_residual_tol
    msg = sprintf(['Density Newton did not converge: ', ...
        'residual=%.3e.'],info.residual);
    if par.abort_on_newton_failure
        error(msg);
    else
        warning('%s',msg);
    end
end
end

%% ------------------------------------------------------------------------
function R = density_residual_2d(n,nold,Ax,Ay,Mmob,par,h)
p = n.^par.gamma;

[Ny,Nx] = size(n);
Fx = zeros(Ny,Nx+1);
Fy = zeros(Ny+1,Nx);
Fx(:,2:Nx) = Ax(:,2:Nx).*(p(:,2:Nx)-p(:,1:Nx-1))/h;
Fy(2:Ny,:) = Ay(2:Ny,:).*(p(2:Ny,:)-p(1:Ny-1,:))/h;

divF = divergence_flux_2d(Fx,Fy,h);
R = n-nold-par.dt*divF+par.dt*Mmob.*(p-par.pH);
end

%% ------------------------------------------------------------------------
function [R,Jmat] = density_residual_jacobian_2d( ...
    n,nold,Ax,Ay,Mmob,par,h)
[Ny,Nx] = size(n);
Ntot = Ny*Nx;
R = density_residual_2d(n,nold,Ax,Ay,Mmob,par,h);

D = par.gamma*n.^(par.gamma-1);
ind = reshape(1:Ntot,Ny,Nx);
coef = par.dt/(h*h);

main = 1+par.dt*Mmob.*D+coef.*D.*( ...
    Ax(:,1:Nx)+Ax(:,2:Nx+1)+Ay(1:Ny,:)+Ay(2:Ny+1,:));

rows = ind(:);
cols = ind(:);
vals = main(:);

if Nx>1
    r = ind(:,2:Nx);
    c = ind(:,1:Nx-1);
    v = -coef*Ax(:,2:Nx).*D(:,1:Nx-1);
    rows = [rows;r(:)]; %#ok<AGROW>
    cols = [cols;c(:)]; %#ok<AGROW>
    vals = [vals;v(:)]; %#ok<AGROW>

    r = ind(:,1:Nx-1);
    c = ind(:,2:Nx);
    v = -coef*Ax(:,2:Nx).*D(:,2:Nx);
    rows = [rows;r(:)]; %#ok<AGROW>
    cols = [cols;c(:)]; %#ok<AGROW>
    vals = [vals;v(:)]; %#ok<AGROW>
end

if Ny>1
    r = ind(2:Ny,:);
    c = ind(1:Ny-1,:);
    v = -coef*Ay(2:Ny,:).*D(1:Ny-1,:);
    rows = [rows;r(:)]; %#ok<AGROW>
    cols = [cols;c(:)]; %#ok<AGROW>
    vals = [vals;v(:)]; %#ok<AGROW>

    r = ind(1:Ny-1,:);
    c = ind(2:Ny,:);
    v = -coef*Ay(2:Ny,:).*D(2:Ny,:);
    rows = [rows;r(:)]; %#ok<AGROW>
    cols = [cols;c(:)]; %#ok<AGROW>
    vals = [vals;v(:)]; %#ok<AGROW>
end

Jmat = sparse(rows,cols,vals,Ntot,Ntot);
end

%% ------------------------------------------------------------------------
function [Ax,Ay] = density_edge_mobilities_2d(nold)
[Ny,Nx] = size(nold);
Ax = zeros(Ny,Nx+1);
Ay = zeros(Ny+1,Nx);
if Nx>1
    Ax(:,2:Nx) = 0.5*(nold(:,1:Nx-1)+nold(:,2:Nx));
end
if Ny>1
    Ay(2:Ny,:) = 0.5*(nold(1:Ny-1,:)+nold(2:Ny,:));
end
end

%% ------------------------------------------------------------------------
function [Jx,Jy] = physical_mass_flux_2d(nold,pnew,h)
%PHYSICAL_MASS_FLUX_2D J=-n^k_edge grad_h p^{k+1}.
[Ny,Nx] = size(nold);
[Ax,Ay] = density_edge_mobilities_2d(nold);
Jx = zeros(Ny,Nx+1);
Jy = zeros(Ny+1,Nx);
if Nx>1
    Jx(:,2:Nx) = -Ax(:,2:Nx).* ...
        (pnew(:,2:Nx)-pnew(:,1:Nx-1))/h;
end
if Ny>1
    Jy(2:Ny,:) = -Ay(2:Ny,:).* ...
        (pnew(2:Ny,:)-pnew(1:Ny-1,:))/h;
end
end

%% ------------------------------------------------------------------------
function divJ = divergence_flux_2d(Jx,Jy,h)
[Ny,Nxp1] = size(Jx);
Nx = Nxp1-1;
if ~isequal(size(Jy),[Ny+1,Nx])
    error('Incompatible x- and y-flux dimensions.');
end
divJ = (Jx(:,2:Nx+1)-Jx(:,1:Nx))/h+ ...
    (Jy(2:Ny+1,:)-Jy(1:Ny,:))/h;
end

%% ------------------------------------------------------------------------
function Mmob = reaction_mobility_species_2d( ...
    n_species,growth_rates)
[Ny,Nx,m] = size(n_species);
if numel(growth_rates)~=m
    error('growth_rates has incompatible length.');
end
Mmob = zeros(Ny,Nx);
for s = 1:m
    Mmob = Mmob+growth_rates(s)*n_species(:,:,s);
end
end

%% ------------------------------------------------------------------------
function [nhat_species,info] = ...
    local_reaction_allocation_diagonal_2d( ...
    n_species_old,nhat,pnew,par)
%LOCAL_REACTION_ALLOCATION_DIAGONAL_2D
% Allocate nhat according to the normalized exact frozen-pressure diagonal
% reaction flow.

[Ny,Nx,m] = size(n_species_old);
nold = sum(n_species_old,3);
y = zeros(Ny,Nx,m);
for s = 1:m
    Hs = par.growth_rates(s)*(par.pH-pnew);
    y(:,:,s) = exp(par.dt*Hs).*n_species_old(:,:,s);
end

den = sum(y,3);
old_active = nold>0;
if any(den(old_active)<=0) || any(~isfinite(den(old_active)))
    error(['Local reaction normalization has a nonpositive or ', ...
        'nonfinite denominator.']);
end

fractions = zeros(Ny,Nx,m);
for s = 1:m
    temp = zeros(Ny,Nx);
    ys = y(:,:,s);
    temp(old_active) = ys(old_active)./den(old_active);
    fractions(:,:,s) = temp;
end

vals = reshape(fractions,Ny*Nx,m);
active_vec = old_active(:);
active_vals = vals(active_vec,:);
if isempty(active_vals)
    overshoot = 0;
else
    overshoot = max([0;-active_vals(:);active_vals(:)-1]);
end
if overshoot>par.remap_bound_tol
    error('Local reaction fractions violate the simplex: %.3e.', ...
        overshoot);
end

fractions = normalize_simplex_on_mask_2d( ...
    fractions,old_active,par.roundoff_tol);

nhat_species = zeros(Ny,Nx,m);
for s = 1:m
    nhat_species(:,:,s) = nhat.*fractions(:,:,s);
end

allocation_array = abs(sum(nhat_species,3)-nhat);
allocation_error = max(allocation_array(:));
info.allocation_error = allocation_error;
info.overshoot = overshoot;
end

%% ------------------------------------------------------------------------
function [n_species,theta,info] = conservative_species_remap_2d( ...
    nnew,nhat_species,Jx,Jy,pnew,dt,h,par)
%CONSERVATIVE_SPECIES_REMAP_2D
% Pressure-ordered solve of the implicit conservative upwind system
%
%   n_c^{k+1} theta_c^{k+1}
%   +dt div_h(J^{k+1} theta_up^{k+1})_c
%   = nhat_c^{k+1}.
%
% Every nonzero total-mass flux points from larger to smaller p^{k+1};
% sorting cells by decreasing pressure therefore yields a triangular sweep.

[Ny,Nx] = size(nnew);
[Ny2,Nx2,m] = size(nhat_species);
if Ny2~=Ny || Nx2~=Nx
    error('nhat_species has incompatible spatial dimensions.');
end
if ~isequal(size(Jx),[Ny,Nx+1]) || ...
        ~isequal(size(Jy),[Ny+1,Nx])
    error('The total-mass flux has incompatible dimensions.');
end
if ~isequal(size(pnew),[Ny,Nx])
    error('pnew has incompatible dimensions.');
end

lambda = dt/h;
Jx_plus = max(Jx,0);
Jx_minus = min(Jx,0);
Jy_plus = max(Jy,0);
Jy_minus = min(Jy,0);

% Diagonal mass: mass remaining in the cell plus all outgoing mass.
d = nnew+lambda*( ...
    Jx_plus(:,2:Nx+1)-Jx_minus(:,1:Nx)+ ...
    Jy_plus(2:Ny+1,:)-Jy_minus(1:Ny,:));

mass_scale = max([1;abs(nnew(:));abs(nhat_species(:)); ...
    lambda*abs(Jx(:));lambda*abs(Jy(:))]);
mass_tol = max(par.remap_mass_tol*mass_scale, ...
    100*eps(mass_scale)*mass_scale);
negative_tol = max(par.roundoff_tol*mass_scale, ...
    100*eps(mass_scale)*mass_scale);

if any(d(:)<-negative_tol)
    error('Remapping diagonal is negative: min(d)=%.3e.',min(d(:)));
end
d(d<0) = 0;
active = d>mass_tol;

Ntot = Ny*Nx;
nhat_mat = reshape(nhat_species,Ntot,m);
theta_mat = zeros(Ntot,m);
processed = false(Ntot,1);
[~,order] = sort(pnew(:),'descend');

max_mass_compatibility_defect = 0;
zero_mass_fallback_count = 0;
max_raw_negative = 0;

for q = 1:Ntot
    idx = order(q);
    if ~active(idx)
        continue;
    end

    iy = mod(idx-1,Ny)+1;
    ix = floor((idx-1)/Ny)+1;
    rhs = nhat_mat(idx,:);

    % Incoming flux from the left cell.
    if ix>1 && lambda*Jx(iy,ix)>mass_tol
        donor = idx-Ny;
        if ~active(donor) || ~processed(donor)
            error('Pressure ordering failed for left inflow into cell %d.',idx);
        end
        rhs = rhs+lambda*Jx(iy,ix)*theta_mat(donor,:);
    end

    % Incoming flux from the right cell.
    if ix<Nx && -lambda*Jx(iy,ix+1)>mass_tol
        donor = idx+Ny;
        if ~active(donor) || ~processed(donor)
            error('Pressure ordering failed for right inflow into cell %d.',idx);
        end
        rhs = rhs-lambda*Jx(iy,ix+1)*theta_mat(donor,:);
    end

    % Incoming flux from the lower cell.
    if iy>1 && lambda*Jy(iy,ix)>mass_tol
        donor = idx-1;
        if ~active(donor) || ~processed(donor)
            error('Pressure ordering failed for lower inflow into cell %d.',idx);
        end
        rhs = rhs+lambda*Jy(iy,ix)*theta_mat(donor,:);
    end

    % Incoming flux from the upper cell.
    if iy<Ny && -lambda*Jy(iy+1,ix)>mass_tol
        donor = idx+1;
        if ~active(donor) || ~processed(donor)
            error('Pressure ordering failed for upper inflow into cell %d.',idx);
        end
        rhs = rhs-lambda*Jy(iy+1,ix)*theta_mat(donor,:);
    end

    local_scale = max([mass_tol,abs(rhs),d(idx)]);
    local_negative_tol = max(par.roundoff_tol*local_scale, ...
        100*eps(local_scale)*local_scale);
    max_raw_negative = max(max_raw_negative,max(0,-min(rhs)));

    if any(rhs<-local_negative_tol)
        error(['Remapping right-hand side is genuinely negative at ', ...
            'cell %d: min=%.3e.'],idx,min(rhs));
    end
    rhs(rhs<0) = 0;
    rhs_mass = sum(rhs);
    max_mass_compatibility_defect = max( ...
        max_mass_compatibility_defect,abs(d(idx)-rhs_mass));

    if rhs_mass>mass_tol
        theta_mat(idx,:) = rhs/rhs_mass;
    else
        local_mass = sum(nhat_mat(idx,:));
        if local_mass>mass_tol
            theta_mat(idx,:) = nhat_mat(idx,:)/local_mass;
        else
            theta_mat(idx,:) = ones(1,m)/m;
        end
        zero_mass_fallback_count = zero_mass_fallback_count+1;
    end

    theta_mat(idx,theta_mat(idx,:)<0) = 0;
    row_sum = sum(theta_mat(idx,:));
    if ~(isfinite(row_sum) && row_sum>0)
        error('Invalid remapping composition at cell %d.',idx);
    end
    theta_mat(idx,:) = theta_mat(idx,:)/row_sum;
    processed(idx) = true;
end

if any(active(:) & ~processed)
    error('Not every active cell was processed in the remapping sweep.');
end

theta_mat(~active(:),:) = 0;
theta = reshape(theta_mat,Ny,Nx,m);

active_vals = theta_mat(active(:),:);
if isempty(active_vals)
    simplex_error = 0;
    overshoot = 0;
    theta_min = 0;
    theta_max = 0;
else
    simplex_error = max(abs(sum(active_vals,2)-1));
    theta_min = min(active_vals(:));
    theta_max = max(active_vals(:));
    overshoot = max([0,-theta_min,theta_max-1]);
end

if overshoot>par.remap_bound_tol || ...
        simplex_error>par.remap_bound_tol
    error(['Normalized remapping composition violates the simplex: ', ...
        'overshoot=%.3e, sum error=%.3e.'], ...
        overshoot,simplex_error);
end

% Reconstruct species densities.  The last component is formed by
% subtraction to preserve sum_s n_s=nnew up to machine arithmetic.
n_species = zeros(Ny,Nx,m);
if m==1
    n_species(:,:,1) = nnew;
else
    for s = 1:m-1
        n_species(:,:,s) = nnew.*theta(:,:,s);
    end
    n_species(:,:,m) = nnew-sum(n_species(:,:,1:m-1),3);
end

species_scale = max(1,max(abs(nnew(:))));
species_negative_tol = par.roundoff_tol*species_scale;
if any(n_species(:)<-species_negative_tol)
    error('Reconstructed species density is negative: min=%.3e.', ...
        min(n_species(:)));
end
n_species(n_species<0) = 0;

% Species fluxes and residual of the implicit conservative system.
Jx_species = zeros(Ny,Nx+1,m);
Jy_species = zeros(Ny+1,Nx,m);
for s = 1:m
    if Nx>1
        th_left = theta(:,1:Nx-1,s);
        th_right = theta(:,2:Nx,s);
        Jint = Jx(:,2:Nx);
        Jx_species(:,2:Nx,s) = ...
            max(Jint,0).*th_left+min(Jint,0).*th_right;
    end
    if Ny>1
        th_lower = theta(1:Ny-1,:,s);
        th_upper = theta(2:Ny,:,s);
        Jint = Jy(2:Ny,:);
        Jy_species(2:Ny,:,s) = ...
            max(Jint,0).*th_lower+min(Jint,0).*th_upper;
    end
end

Rspecies = zeros(Ny,Nx,m);
for s = 1:m
    divJs = divergence_flux_2d( ...
        Jx_species(:,:,s),Jy_species(:,:,s),h);
    Rspecies(:,:,s) = n_species(:,:,s)+dt*divJs- ...
        nhat_species(:,:,s);
end

reconstruction_array = abs(sum(n_species,3)-nnew);
reconstruction_error = max(reconstruction_array(:));
species_residual = max(abs(Rspecies(:)));
residual_scale = 1+max(abs(nhat_species(:)));

if reconstruction_error>par.reconstruction_tol
    error('Species reconstruction error is %.3e.', ...
        reconstruction_error);
end
if species_residual>par.remap_residual_tol*residual_scale
    msg = sprintf(['Species remapping residual is %.3e; ', ...
        'scalar/species mass-compatibility defect is %.3e.'], ...
        species_residual,max_mass_compatibility_defect);
    if par.abort_on_remap_failure
        error(msg);
    else
        warning('%s',msg);
    end
end

info.active_cell_count = nnz(active);
info.raw_min_theta = theta_min;
info.raw_max_theta = theta_max;
info.overshoot = overshoot;
info.simplex_error = simplex_error;
info.reconstruction_error = reconstruction_error;
info.species_residual = species_residual;
info.min_diagonal = min_positive_2d(d(active));
info.mass_compatibility_defect = max_mass_compatibility_defect;
info.zero_mass_fallback_count = zero_mass_fallback_count;
info.max_raw_negative_rhs = max_raw_negative;
info.mass_tolerance = mass_tol;
end

%% ------------------------------------------------------------------------
function [Dtrans,Dreact,defect] = density_energy_diagnostics_2d( ...
    nold,pnew,Mmob,Eold,Enew,par,h)
[Ny,Nx] = size(nold);
[Ax,Ay] = density_edge_mobilities_2d(nold);

Dtrans = 0;
if Nx>1
    gradx = (pnew(:,2:Nx)-pnew(:,1:Nx-1))/h;
    temp = Ax(:,2:Nx).*gradx.^2;
    Dtrans = Dtrans+par.dt*h^2*sum(temp(:));
end
if Ny>1
    grady = (pnew(2:Ny,:)-pnew(1:Ny-1,:))/h;
    temp = Ay(2:Ny,:).*grady.^2;
    Dtrans = Dtrans+par.dt*h^2*sum(temp(:));
end
temp = Mmob.*(pnew-par.pH).^2;
Dreact = par.dt*h^2*sum(temp(:));
defect = Enew+Dtrans+Dreact-Eold;
end

%% ------------------------------------------------------------------------
function [diag_out,species_masses] = state_diagnostics_2d( ...
    n_species,n,p,h,par)
reconstruction = sum(n_species,3);
phi = species_fractions_2d(n_species,n,par.occupied_tol);
occupied = n>par.occupied_tol;

if any(occupied(:))
    phi_mat = reshape(phi,[],size(phi,3));
    occupied_values = phi_mat(occupied(:),:);
    simplex_error = max(abs(sum(occupied_values,2)-1));
else
    simplex_error = 0;
end

diag_out.min_species = min(n_species(:));
diag_out.min_n = min(n(:));
diag_out.max_n = max(n(:));
diag_out.max_p = max(p(:));
diag_out.reconstruction_error = max(abs(reconstruction(:)-n(:)));
diag_out.simplex_error = simplex_error;

m = size(n_species,3);
species_masses = zeros(1,m);
for s = 1:m
    temp = n_species(:,:,s);
    species_masses(s) = h^2*sum(temp(:));
end
end

%% ------------------------------------------------------------------------
function phi = species_fractions_2d(n_species,n,occupied_tol)
[Ny,Nx,m] = size(n_species);
phi = NaN(Ny,Nx,m);
occupied = n>occupied_tol;
for s = 1:m
    temp = NaN(Ny,Nx);
    ns = n_species(:,:,s);
    temp(occupied) = ns(occupied)./n(occupied);
    phi(:,:,s) = temp;
end
end

%% ------------------------------------------------------------------------
function phi = normalize_simplex_array_2d(phi,roundoff_tol)
[Ny,Nx,m] = size(phi);
phi(phi<0 & phi>=-roundoff_tol) = 0;
if any(phi(:)<0)
    error('Simplex array contains genuinely negative entries.');
end
S = sum(phi,3);
if any(~isfinite(S(:)) | S(:)<=0)
    error('Simplex normalization has a nonpositive or nonfinite sum.');
end
for s = 1:m
    phi(:,:,s) = phi(:,:,s)./S;
end
end

%% ------------------------------------------------------------------------
function phi = normalize_simplex_on_mask_2d(phi,mask,roundoff_tol)
[Ny,Nx,m] = size(phi);
phi(phi<0 & phi>=-roundoff_tol) = 0;
if any(phi(:)<0)
    error('Simplex entries are genuinely negative.');
end
S = sum(phi,3);
if any(S(mask)<=0) || any(~isfinite(S(mask)))
    error('Active simplex normalization has an invalid sum.');
end
for s = 1:m
    temp = phi(:,:,s);
    temp(mask) = temp(mask)./S(mask);
    temp(~mask) = 0;
    phi(:,:,s) = temp;
end
end

%% ------------------------------------------------------------------------
function y = roundoff_project_interval_2d( ...
    y,lo,hi,roundoff_tol,label)
low_round = y<lo & y>=lo-roundoff_tol;
high_round = y>hi & y<=hi+roundoff_tol;
y(low_round) = lo;
y(high_round) = hi;

if any(y(:)<lo | y(:)>hi)
    warning('%s exceeds [%g,%g]: min=%.3e, max=%.3e.', ...
        label,lo,hi,min(y(:)),max(y(:)));
end
end

%% ------------------------------------------------------------------------
function E = modified_energy_2d(n,h,gamma,pH)
density = n.^(gamma+1)/(gamma+1)-pH*n;
E = h^2*sum(density(:));
end

%% ------------------------------------------------------------------------
function value = min_positive_2d(x)
pos = x(x>0);
if isempty(pos)
    value = NaN;
else
    value = min(pos);
end
end

%% ------------------------------------------------------------------------
function snaps = store_snapshot_2d( ...
    snaps,sid,t,n,p,n_species)
snaps.t(sid) = t;
snaps.n{sid} = n;
snaps.p{sid} = p;
snaps.n_species{sid} = n_species;
end

%% ------------------------------------------------------------------------
function [theta_vec,r_front] = compute_front_radius_2d( ...
    x,y,n,front_level,n_theta)
theta_vec = linspace(-pi,pi,n_theta);
r_front = NaN(size(theta_vec));

r_max = 0.98*min(max(abs(x)),max(abs(y)));
nr = 1200;
r_line = linspace(0,r_max,nr);

for q = 1:numel(theta_vec)
    th = theta_vec(q);
    xq = r_line*cos(th);
    yq = r_line*sin(th);
    vals = interp2(x,y,n,xq,yq,'linear',0);
    idx = find(vals>=front_level,1,'last');
    if isempty(idx)
        continue;
    elseif idx>=numel(r_line)
        r_front(q) = r_line(end);
    else
        r1 = r_line(idx);
        r2 = r_line(idx+1);
        v1 = vals(idx);
        v2 = vals(idx+1);
        if abs(v2-v1)<1.0e-14
            r_front(q) = r1;
        else
            r_front(q) = r1+(front_level-v1)* ...
                (r2-r1)/(v2-v1);
        end
    end
end
end

%% ------------------------------------------------------------------------
