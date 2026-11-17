function run_example2_transition_1d(outdir)
%RUN_EXAMPLE2_TRANSITION_1D
% Example 2: pressure-induced phenotype switching with a common total
% growth law.
%
% Numerical scheme:
%   1. scalar Onsager density-pressure step with frozen mobility;
%   2. local normalized exact frozen-pressure reaction allocation;
%   3. implicit conservative upwind species remapping driven by exactly
%      the same physical mass flux as the total-density step.
%
% The fully discrete state is (n_1,n_2).  The composition is reconstructed
% only in occupied cells and is not transported as an independent all-grid
% variable.  The remapping preserves species nonnegativity and
% n_1+n_2=n without an additional species CFL condition.
%
% This file is self-contained and compatible with older MATLAB versions.

clc; close all;
if nargin < 1 || isempty(outdir), outdir=fullfile(tempdir,'two_species_rerun','switching'); end

%% Parameters
par=model.switching_parameters;
fprintf('Code version: %s\n',par.code_version);

c_star = par.c_star;

%% Output folder
% outdir is supplied by the caller; numerical parameters are unchanged.
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% Grid and time
N = par.N;
h = (par.b-par.a)/N;
x = par.a + ((1:N)'-0.5)*h;
[~,center_idx] = min(abs(x));

nt = round(par.T/par.dt);
par.dt = par.T/nt;
time_grid = (0:nt)'*par.dt;

snapshot_steps = round(par.snapshot_times/par.dt);
snapshot_steps = max(0,min(nt,snapshot_steps));
[snapshot_steps,ia] = unique(snapshot_steps,'stable');
snapshot_times = snapshot_steps(:)*par.dt;
par.snapshot_times = par.snapshot_times(ia);
nsnapshot = numel(snapshot_steps);

assert(~isempty(snapshot_steps) && snapshot_steps(1) == 0, ...
    'snapshot_steps must include the initial time.');
assert(snapshot_steps(end) == nt, ...
    'snapshot_steps must include the final time.');

%% Initial data
% Same total-density profile as in the structure-preservation and
% non-closure tests:
%   n^0(x)=0.49[tanh((x+2)/0.05)-tanh((x-2)/0.05)].
[n0,phi1_0,phi2_0,n1_0,n2_0]=model.switching_initial(x,par);

fprintf('Running Example 2 with the conservative species remapping.\n');
result = run_single_case(n1_0,n2_0,x,h,nt,snapshot_steps,center_idx,par);

%% Backward-compatible aliases
% The old code used c=phi_1.  These aliases allow older plotting scripts to
% read the new data file, while the actual discrete state remains (n_1,n_2).
t_hist = snapshot_times;
n_hist = result.n_hist;
p_hist = result.p_hist;
phi1_hist = result.phi1_hist;
phi2_hist = result.phi2_hist;
c_hist = phi1_hist;
n1_hist = result.n1_hist;
n2_hist = result.n2_hist;

mass_hist = result.mass_hist;
mass1_hist = result.mass1_hist;
mass2_hist = result.mass2_hist;
fraction1_hist = result.fraction1_hist;
fraction2_hist = result.fraction2_hist;
mean_c_hist = fraction1_hist;
center_phi1_hist = result.center_phi1_hist;
center_phi2_hist = result.center_phi2_hist;
core_phi1_hist = result.core_phi1_hist;
core_phi2_hist = result.core_phi2_hist;
front_phi1_hist = result.front_phi1_hist;
front_phi2_hist = result.front_phi2_hist;
center_c_hist = center_phi1_hist;
core_c_hist = core_phi1_hist;
front_c_hist = front_phi1_hist;
energy_hist = result.energy_hist;
newton_it_hist = result.newton_it_hist;
newton_res_hist = result.newton_res_hist;
remap_overshoot_count = result.remap_overshoot_count;
remap_overshoot_max = result.remap_overshoot_max;

%% Save data
save(fullfile(outdir,'example2_data_conservative_remap.mat'), ...
    'x','t_hist','snapshot_times','time_grid', ...
    'n0','phi1_0','phi2_0','n1_0','n2_0', ...
    'n_hist','p_hist','phi1_hist','phi2_hist','c_hist', ...
    'n1_hist','n2_hist', ...
    'mass_hist','mass1_hist','mass2_hist', ...
    'fraction1_hist','fraction2_hist','mean_c_hist', ...
    'center_phi1_hist','center_phi2_hist', ...
    'core_phi1_hist','core_phi2_hist', ...
    'front_phi1_hist','front_phi2_hist', ...
    'center_c_hist','core_c_hist','front_c_hist', ...
    'center_idx','c_star','energy_hist', ...
    'newton_it_hist','newton_res_hist', ...
    'remap_overshoot_count','remap_overshoot_max', ...
    'result','par','-v7.3');
print_case_summary(result,par);

fprintf('\nExample 2 finished.\n');
fprintf('Output folder: %s\n',outdir);

end

%% ------------------------------------------------------------------------
function result = run_single_case( ...
    n1_0,n2_0,x,h,nt,snapshot_steps,center_idx,par)
N = length(n1_0);
nsnapshot = numel(snapshot_steps);

n1 = max(n1_0,0);
n2 = max(n2_0,0);
n = n1+n2;
p = n.^par.gamma;
phi1 = species_fraction_1(n1,n2,par.occupied_tol);
phi2 = species_fraction_2(n1,n2,par.occupied_tol);

n_hist = zeros(N,nsnapshot);
p_hist = zeros(N,nsnapshot);
phi1_hist = NaN(N,nsnapshot);
phi2_hist = NaN(N,nsnapshot);
n1_hist = zeros(N,nsnapshot);
n2_hist = zeros(N,nsnapshot);

mass_hist = zeros(nt+1,1);
mass1_hist = zeros(nt+1,1);
mass2_hist = zeros(nt+1,1);
fraction1_hist = zeros(nt+1,1);
fraction2_hist = zeros(nt+1,1);
center_phi1_hist = NaN(nt+1,1);
center_phi2_hist = NaN(nt+1,1);
core_phi1_hist = NaN(nt+1,1);
core_phi2_hist = NaN(nt+1,1);
front_phi1_hist = NaN(nt+1,1);
front_phi2_hist = NaN(nt+1,1);
energy_hist = zeros(nt+1,1);

newton_it_hist = zeros(nt,1);
newton_res_hist = zeros(nt,1);
transport_diss_hist = zeros(nt,1);
reaction_diss_hist = zeros(nt,1);
energy_ineq_defect_hist = zeros(nt,1);
source_balance_defect_hist = zeros(nt,1);
reaction_allocation_err_hist = zeros(nt,1);
reaction_fraction_overshoot_hist = zeros(nt,1);
remap_simplex_err_hist = zeros(nt,1);
remap_reconstruction_err_hist = zeros(nt,1);
remap_residual_hist = zeros(nt,1);
active_cell_count_hist = zeros(nt,1);

remap_overshoot_count = 0;
remap_overshoot_max = 0.0;

snapshot_id = 1;
if snapshot_id <= nsnapshot && snapshot_steps(snapshot_id) == 0
    [n_hist,p_hist,phi1_hist,phi2_hist,n1_hist,n2_hist] = ...
        store_snapshot(snapshot_id,n,p,phi1,phi2,n1,n2, ...
        n_hist,p_hist,phi1_hist,phi2_hist,n1_hist,n2_hist);
    snapshot_id = snapshot_id+1;
end

[mass_hist(1),mass1_hist(1),mass2_hist(1), ...
    fraction1_hist(1),fraction2_hist(1)] = ...
    compute_masses(n,n1,n2,h);
center_phi1_hist(1) = phi1(center_idx);
center_phi2_hist(1) = phi2(center_idx);
[core_phi1_hist(1),front_phi1_hist(1)] = ...
    compute_core_front_averages(n,p,phi1,par);
[core_phi2_hist(1),front_phi2_hist(1)] = ...
    compute_core_front_averages(n,p,phi2,par);
energy_hist(1) = energy_h(n,par,h);

theta = zeros(N,2);

for k = 1:nt
    nold = n;
    n1old = n1;
    n2old = n2;
    Eold = energy_hist(k);

    old_reconstruction_err = max(abs(n1old+n2old-nold));
    if old_reconstruction_err > par.reconstruction_tol
        error('Reconstruction error before step %d is %.3e.', ...
            k,old_reconstruction_err);
    end

    % Since H_1=H_2=a_H(p_H-p), the scalar reaction mobility is a_H*n^k.
    Mmob = reaction_mobility_species(n1old,n2old,par);

    % Step 1: scalar Onsager density-pressure update.
    [nnew,pnew,newton_info] = solve_density_step_n(nold,Mmob,par,h);
    newton_it_hist(k) = newton_info.iter;
    newton_res_hist(k) = newton_info.residual;

    Enew = energy_h(nnew,par,h);
    [Dtrans,Dreact,Edefect] = density_energy_diagnostics( ...
        nold,pnew,Mmob,Eold,Enew,par,h);
    transport_diss_hist(k) = Dtrans;
    reaction_diss_hist(k) = Dreact;
    energy_ineq_defect_hist(k) = Edefect;

    if Edefect > par.energy_ineq_tol
        warning('Discrete energy inequality defect at step %d is %.3e.', ...
            k,Edefect);
    end

    Jmass = physical_mass_flux(nold,pnew,h);

    % The source and conservative forms coincide for an exact nonlinear
    % solve.  Their difference is retained as a Newton diagnostic.
    S = Mmob.*(par.pH-pnew);
    source_tol = par.roundoff_tol*(1+max(abs(S)));
    if any(S < -source_tol)
        error('Negative density source at step %d: min(S)=%.3e.',k,min(S));
    end
    S(S < 0) = 0;

    nhat_source = nold+par.dt*S;
    nhat_balance = nnew+(par.dt/h)*(Jmass(2:end)-Jmass(1:end-1));
    source_balance_defect = max(abs(nhat_balance-nhat_source));
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

    nhat_scale = max(1.0,max(abs(nhat)));
    nhat_tol = max( ...
        par.roundoff_tol*nhat_scale, ...
        100*eps(nhat_scale)*nhat_scale);
    small_negative = nhat < 0 & nhat >= -nhat_tol;
    nhat(small_negative) = 0;
    if any(nhat < -nhat_tol)
        error('Physical source-updated mass is negative at step %d: min=%.3e.', ...
            k,min(nhat));
    end

    % Step 2: local normalized exact frozen-pressure reaction allocation.
    [nhat_species,reaction_info] = local_reaction_allocation_example2( ...
        n1old,n2old,nhat,pnew,par);
    reaction_allocation_err_hist(k) = reaction_info.allocation_error;
    reaction_fraction_overshoot_hist(k) = reaction_info.overshoot;

    % Step 3: implicit conservative upwind species remapping driven by J.
    [n_species,theta,remap_info] = conservative_species_remap_1d( ...
        nnew,nhat_species,Jmass,pnew,par.dt,h,par);

    n1 = n_species(:,1);
    n2 = n_species(:,2);
    n = nnew;
    p = pnew;
    phi1 = species_fraction_1(n1,n2,par.occupied_tol);
    phi2 = species_fraction_2(n1,n2,par.occupied_tol);

    remap_simplex_err_hist(k) = remap_info.simplex_error;
    remap_reconstruction_err_hist(k) = remap_info.reconstruction_error;
    remap_residual_hist(k) = remap_info.species_residual;
    active_cell_count_hist(k) = remap_info.active_cell_count;

    if remap_info.overshoot > par.remap_bound_tol
        remap_overshoot_count = remap_overshoot_count+1;
    end
    remap_overshoot_max = max(remap_overshoot_max,remap_info.overshoot);

    [mass_hist(k+1),mass1_hist(k+1),mass2_hist(k+1), ...
        fraction1_hist(k+1),fraction2_hist(k+1)] = ...
        compute_masses(n,n1,n2,h);
    center_phi1_hist(k+1) = phi1(center_idx);
    center_phi2_hist(k+1) = phi2(center_idx);
    [core_phi1_hist(k+1),front_phi1_hist(k+1)] = ...
        compute_core_front_averages(n,p,phi1,par);
    [core_phi2_hist(k+1),front_phi2_hist(k+1)] = ...
        compute_core_front_averages(n,p,phi2,par);
    energy_hist(k+1) = Enew;

    if snapshot_id <= nsnapshot && k == snapshot_steps(snapshot_id)
        [n_hist,p_hist,phi1_hist,phi2_hist,n1_hist,n2_hist] = ...
            store_snapshot(snapshot_id,n,p,phi1,phi2,n1,n2, ...
            n_hist,p_hist,phi1_hist,phi2_hist,n1_hist,n2_hist);
        snapshot_id = snapshot_id+1;
    end

    if mod(k,par.progress_every) == 0 || k == nt
        fprintf(['step %5d/%5d, t=%.4f, Newton it=%2d, res=%.3e, ', ...
            'energy defect=%.3e, remap res=%.3e, recon=%.3e\n'], ...
            k,nt,k*par.dt,newton_info.iter,newton_info.residual, ...
            Edefect,remap_info.species_residual, ...
            remap_info.reconstruction_error);
    end
end

if snapshot_id ~= nsnapshot+1
    error('Not every requested Example 2 snapshot was saved.');
end

result.n_hist = n_hist;
result.p_hist = p_hist;
result.phi1_hist = phi1_hist;
result.phi2_hist = phi2_hist;
result.n1_hist = n1_hist;
result.n2_hist = n2_hist;
result.mass_hist = mass_hist;
result.mass1_hist = mass1_hist;
result.mass2_hist = mass2_hist;
result.fraction1_hist = fraction1_hist;
result.fraction2_hist = fraction2_hist;
result.center_phi1_hist = center_phi1_hist;
result.center_phi2_hist = center_phi2_hist;
result.core_phi1_hist = core_phi1_hist;
result.core_phi2_hist = core_phi2_hist;
result.front_phi1_hist = front_phi1_hist;
result.front_phi2_hist = front_phi2_hist;
result.energy_hist = energy_hist;
result.newton_it_hist = newton_it_hist;
result.newton_res_hist = newton_res_hist;
result.transport_diss_hist = transport_diss_hist;
result.reaction_diss_hist = reaction_diss_hist;
result.energy_ineq_defect_hist = energy_ineq_defect_hist;
result.source_balance_defect_hist = source_balance_defect_hist;
result.reaction_allocation_err_hist = reaction_allocation_err_hist;
result.reaction_fraction_overshoot_hist = reaction_fraction_overshoot_hist;
result.remap_simplex_err_hist = remap_simplex_err_hist;
result.remap_reconstruction_err_hist = remap_reconstruction_err_hist;
result.remap_residual_hist = remap_residual_hist;
result.active_cell_count_hist = active_cell_count_hist;
result.remap_overshoot_count = remap_overshoot_count;
result.remap_overshoot_max = remap_overshoot_max;
result.final_n = n;
result.final_p = p;
result.final_phi1 = phi1;
result.final_phi2 = phi2;
result.final_n1 = n1;
result.final_n2 = n2;
result.final_theta1 = theta(:,1);
result.final_theta2 = theta(:,2);
end

function [nnew,pnew,info] = solve_density_step_n(nold,Mmob,par,h)
%SOLVE_DENSITY_STEP_N Projected damped Newton solve of the scalar
% density-pressure Onsager step.
%
% The equation is
%   n^{k+1}-n^k
%   = dt*d_h(n^k_edge*delta_h p^{k+1})
%     -dt*M^k*(p^{k+1}-p_H),
%   p^{k+1}=(n^{k+1})^gamma.
%
% This is the same nonlinear scalar solve as in the single-population
% method.  The only change is that M^k is composition dependent.

N = length(nold);
n = min(max(nold,0),par.nH);

% Interface coefficients n_{j+1/2}^k.
edge_n = zeros(N+1,1);
if N > 1
    edge_n(2:N) = 0.5*(nold(1:N-1)+nold(2:N));
end

info.converged = false;
info.iter = 0;
info.residual = Inf;
info.total_line_search_steps = 0;
info.min_accepted_lambda = 1.0;

stop_tol = par.newton_tol*(1+norm(nold,inf));

for it = 1:par.newton_maxit
    [R,Jmat] = density_residual_jac(n,nold,edge_n,Mmob,par,h);
    res_inf = norm(R,inf);

    if res_inf <= stop_tol
        info.converged = true;
        info.iter = it-1;
        info.residual = res_inf;
        break;
    end

    dn = -(Jmat\R);
    if any(~isfinite(dn))
        error('Density Newton produced a non-finite search direction.');
    end

    lambda = 1.0;
    accepted = false;
    res2 = norm(R,2);
    best_n = n;
    best_res2 = res2;
    best_ls = 0;

    for ls = 1:par.line_search_maxit
        % Projected trial step keeps every Newton trial in the invariant box.
        % The lower bound is exactly zero, so genuine vacuum is retained.
        ntrial = n + lambda*dn;
        ntrial = min(max(ntrial,0),par.nH);

        Rtrial = density_residual(ntrial,nold,edge_n,Mmob,par,h);
        trial_res2 = norm(Rtrial,2);

        if isfinite(trial_res2) && trial_res2 < best_res2
            best_res2 = trial_res2;
            best_n = ntrial;
            best_ls = ls;
        end

        if isfinite(trial_res2) && trial_res2 <= (1-par.armijo*lambda)*res2
            accepted = true;
            break;
        end

        lambda = 0.5*lambda;
    end

    info.total_line_search_steps = info.total_line_search_steps + ls;

    if accepted
        n = ntrial;
        info.iter = it;
        info.residual = norm(Rtrial,inf);
        info.min_accepted_lambda = min(info.min_accepted_lambda,lambda);
    elseif best_res2 < res2
        n = best_n;
        info.iter = it;
        info.residual = norm(density_residual(n,nold,edge_n,Mmob,par,h),inf);
        if best_ls > 0
            info.min_accepted_lambda = min(info.min_accepted_lambda,2^(-(best_ls-1)));
        end
    else
        info.iter = it;
        info.residual = res_inf;
        warning('Density Newton line search failed at iteration %d; residual = %.3e.', ...
            it,res_inf);
        break;
    end
end

nnew = roundoff_project_interval(n,0,par.nH, ...
    par.roundoff_tol,'density n');
pnew = nnew.^par.gamma;

Rfinal = density_residual(nnew,nold,edge_n,Mmob,par,h);
info.residual = norm(Rfinal,inf);
if info.residual <= stop_tol
    info.converged = true;
end

if ~info.converged && info.residual > par.warn_residual_tol
    msg = sprintf('Density Newton did not converge: residual = %.3e.',info.residual);
    if par.abort_on_newton_failure
        error(msg);
    else
        warning('%s',msg);
    end
end

end

%% ------------------------------------------------------------------------
function R = density_residual(n,nold,edge_n,Mmob,par,h)
N = length(n);
p = n.^par.gamma;

F = zeros(N+1,1);
if N > 1
    F(2:N) = edge_n(2:N).*(p(2:N)-p(1:N-1))/h;
end

dhF = (F(2:N+1)-F(1:N))/h;
R = n - nold - par.dt*dhF + par.dt*Mmob.*(p-par.pH);
end

%% ------------------------------------------------------------------------
function [R,Jmat] = density_residual_jac(n,nold,edge_n,Mmob,par,h)
N = length(n);
R = density_residual(n,nold,edge_n,Mmob,par,h);

dpdn = par.gamma*n.^(par.gamma-1);
coef = par.dt/(h*h);

row = zeros(3*N,1);
col = zeros(3*N,1);
val = zeros(3*N,1);
cnt = 0;

for j = 1:N
    Aleft = edge_n(j);
    Aright = edge_n(j+1);

    diag_entry = 1 ...
        + par.dt*Mmob(j)*dpdn(j) ...
        + coef*(Aleft+Aright)*dpdn(j);

    cnt = cnt + 1;
    row(cnt) = j;
    col(cnt) = j;
    val(cnt) = diag_entry;

    if j > 1
        cnt = cnt + 1;
        row(cnt) = j;
        col(cnt) = j-1;
        val(cnt) = -coef*Aleft*dpdn(j-1);
    end

    if j < N
        cnt = cnt + 1;
        row(cnt) = j;
        col(cnt) = j+1;
        val(cnt) = -coef*Aright*dpdn(j+1);
    end
end

Jmat = sparse(row(1:cnt),col(1:cnt),val(1:cnt),N,N);
end

%% ------------------------------------------------------------------------

%% ------------------------------------------------------------------------
function J = physical_mass_flux(nold,pnew,h)
%PHYSICAL_MASS_FLUX Physical total-density flux J=-n^k grad_h p^{k+1}.
N = length(nold);
F = zeros(N+1,1);
if N > 1
    edge_n = 0.5*(nold(1:N-1)+nold(2:N));
    F(2:N) = edge_n.*(pnew(2:N)-pnew(1:N-1))/h;
end
J = -F;
end

%% ------------------------------------------------------------------------

%% ------------------------------------------------------------------------
function [nhat_species,info] = local_reaction_allocation_example2( ...
    n1old,n2old,nhat,pnew,par)
%LOCAL_REACTION_ALLOCATION_EXAMPLE2
% Allocate the source-updated total mass by the normalized exact
% frozen-pressure reaction flow.
%
% With
%   u = alpha*p,             type 1 -> type 2,
%   v = beta*(p_H-p),        type 2 -> type 1,
% the common total growth a_H(p_H-p) cancels from the normalized reaction
% map.  The type-1 fraction therefore satisfies exactly
%   phi_1^{react} = phi_1^* + (phi_1^old-phi_1^*) exp[-dt(u+v)],
% where phi_1^*=v/(u+v).

nold = n1old+n2old;
old_active = nold > 0;

phi1old = zeros(size(nold));
phi1old(old_active) = n1old(old_active)./nold(old_active);

u = par.alpha*pnew;
v = par.beta*(par.pH-pnew);
rate = u+v;

rate_tol = 100*eps(max(1,max(abs(rate))));
phi1react = phi1old;
positive_rate = old_active & rate > rate_tol;

phi1eq = zeros(size(rate));
phi1eq(positive_rate) = v(positive_rate)./rate(positive_rate);
phi1react(positive_rate) = phi1eq(positive_rate) + ...
    (phi1old(positive_rate)-phi1eq(positive_rate)) ...
    .*exp(-par.dt*rate(positive_rate));

values = phi1react(old_active);
if isempty(values)
    overshoot = 0;
else
    overshoot = max([0;-values;values-1]);
end
if overshoot > par.remap_bound_tol
    error('Local reaction fraction violates [0,1]: overshoot=%.3e.',overshoot);
end
phi1react(old_active) = roundoff_project_interval( ...
    phi1react(old_active),0,1,par.roundoff_tol, ...
    'local reaction fraction');
phi1react(~old_active) = 0;

nhat1 = nhat.*phi1react;
nhat2 = nhat-nhat1;
nhat_species = [nhat1,nhat2];

info.allocation_error = max(abs(sum(nhat_species,2)-nhat));
info.overshoot = overshoot;
end

function [n_species,theta,info] = conservative_species_remap_1d( ...
    nnew,nhat_species,J,pnew,dt,h,par)
%CONSERVATIVE_SPECIES_REMAP_1D
% Robust pressure-ordered solve of the implicit conservative upwind system.
%
% The exact equations imply that the two components of each active-cell
% right-hand side sum to the scalar diagonal mass d_j.  In floating-point
% arithmetic the scalar Newton solve leaves a tiny compatibility defect,
% especially near vacuum.  We therefore normalize the nonnegative species
% right-hand side by its own component sum.  This is algebraically identical
% to division by d_j in exact arithmetic and enforces the simplex constraint
% without changing the total-density solution.

N = length(nnew);
if size(nhat_species,1) ~= N || size(nhat_species,2) ~= 2
    error('The two-species remapping expects an N-by-2 species array.');
end
if length(J) ~= N+1 || length(pnew) ~= N
    error('Incompatible scalar flux or pressure dimensions in remapping.');
end

lambda = dt/h;
Jplus = max(J,0);
Jminus = min(J,0);

d = nnew + lambda*(Jplus(2:N+1)-Jminus(1:N));
mass_scale = max([1; abs(nnew(:)); abs(nhat_species(:)); lambda*abs(J(:))]);
mass_tol = max(par.remap_mass_tol*mass_scale,100*eps(mass_scale)*mass_scale);
negative_tol = max(par.roundoff_tol*mass_scale,100*eps(mass_scale)*mass_scale);

if any(d < -negative_tol)
    error('Remapping diagonal is negative: min(d)=%.3e.',min(d));
end
d(d < 0) = 0;
active = d > mass_tol;

[~,order] = sort(pnew,'descend');
theta = zeros(N,2);
processed = false(N,1);
max_mass_compatibility_defect = 0;
zero_mass_fallback_count = 0;
max_raw_negative = 0;

for q = 1:N
    j = order(q);
    if ~active(j)
        continue;
    end

    rhs = nhat_species(j,:);

    if j > 1 && J(j) > 0
        if ~active(j-1) || ~processed(j-1)
            error('Pressure ordering failed for left inflow into cell %d.',j);
        end
        rhs = rhs + lambda*J(j)*theta(j-1,:);
    end

    if j < N && J(j+1) < 0
        if ~active(j+1) || ~processed(j+1)
            error('Pressure ordering failed for right inflow into cell %d.',j);
        end
        rhs = rhs - lambda*J(j+1)*theta(j+1,:);
    end

    local_scale = max([mass_tol,abs(rhs),d(j)]);
    local_neg_tol = max(par.roundoff_tol*local_scale,100*eps(local_scale)*local_scale);
    max_raw_negative = max(max_raw_negative,max(0,-min(rhs)));
    if any(rhs < -local_neg_tol)
        error('Remapping right-hand side is genuinely negative at cell %d: min = %.3e.', ...
            j,min(rhs));
    end
    rhs(rhs < 0) = 0;
    rhs_mass = sum(rhs);
    max_mass_compatibility_defect = max( ...
        max_mass_compatibility_defect,abs(d(j)-rhs_mass));

    if rhs_mass > mass_tol
        theta(j,:) = rhs/rhs_mass;
    else
        % Numerical vacuum: the carried mass is below the resolution used by
        % the remapping diagnostics.  Its composition has no observable
        % effect on n_i=n*theta.  Use the local old/reaction fraction if
        % available, otherwise a neutral fraction.
        local_mass = sum(nhat_species(j,:));
        if local_mass > mass_tol
            theta(j,:) = nhat_species(j,:)/local_mass;
        else
            theta(j,:) = [0.5,0.5];
        end
        zero_mass_fallback_count = zero_mass_fallback_count + 1;
    end

    % Enforce the simplex locally before this composition is sent downstream.
    theta(j,theta(j,:) < 0) = 0;
    row_sum = sum(theta(j,:));
    if ~(isfinite(row_sum) && row_sum > 0)
        error('Invalid remapping composition at cell %d.',j);
    end
    theta(j,:) = theta(j,:)/row_sum;
    processed(j) = true;
end

if any(active & ~processed)
    error('Not every active cell was processed in the remapping sweep.');
end

theta(~active,:) = 0;
vals = theta(active,:);
if isempty(vals)
    simplex_error = 0;
    overshoot = 0;
    theta_min = 0;
    theta_max = 0;
else
    simplex_error = max(abs(sum(vals,2)-1));
    theta_min = min(vals(:));
    theta_max = max(vals(:));
    overshoot = max([0,-theta_min,theta_max-1]);
end

% This check is applied after local normalization.  A raw compatibility
% defect is recorded separately and is not itself a simplex violation.
if overshoot > par.remap_bound_tol || simplex_error > par.remap_bound_tol
    error(['Normalized remapping composition violates the simplex: ', ...
        'overshoot = %.3e, sum error = %.3e.'],overshoot,simplex_error);
end

n1new = nnew.*theta(:,1);
n2new = nnew-n1new;
n_species = [n1new,n2new];

Jspecies = zeros(N+1,2);
for e = 2:N
    if J(e) > 0
        Jspecies(e,:) = J(e)*theta(e-1,:);
    elseif J(e) < 0
        Jspecies(e,:) = J(e)*theta(e,:);
    end
end
Rspecies = n_species + lambda*(Jspecies(2:N+1,:)-Jspecies(1:N,:)) ...
    - nhat_species;

reconstruction_error = max(abs(sum(n_species,2)-nnew));
species_residual = max(abs(Rspecies(:)));
residual_scale = 1+max(abs(nhat_species(:)));

if reconstruction_error > par.reconstruction_tol
    error('Species reconstruction error is %.3e.',reconstruction_error);
end
if species_residual > par.remap_residual_tol*residual_scale
    msg = sprintf(['Species remapping residual is %.3e; scalar/species ', ...
        'mass-compatibility defect is %.3e.'], ...
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
info.min_diagonal = min_positive(d(active));
info.mass_compatibility_defect = max_mass_compatibility_defect;
info.zero_mass_fallback_count = zero_mass_fallback_count;
info.max_raw_negative_rhs = max_raw_negative;
info.mass_tolerance = mass_tol;
end

%% ------------------------------------------------------------------------
function value = min_positive(x)
pos = x(x>0);
if isempty(pos)
    value = NaN;
else
    value = min(pos);
end
end


%% ------------------------------------------------------------------------
function Mmob = reaction_mobility_species(n1,n2,par)
% Because H_1=H_2=a_H(p_H-p), the frozen scalar mobility is a_H*n^k.
Mmob = par.aH*(n1+n2);
end

%% ------------------------------------------------------------------------
function phi1 = species_fraction_1(n1,n2,occupied_tol)
n = n1+n2;
phi1 = NaN(size(n));
occupied = n > occupied_tol;
phi1(occupied) = n1(occupied)./n(occupied);
end

%% ------------------------------------------------------------------------
function phi2 = species_fraction_2(n1,n2,occupied_tol)
n = n1+n2;
phi2 = NaN(size(n));
occupied = n > occupied_tol;
phi2(occupied) = n2(occupied)./n(occupied);
end

%% ------------------------------------------------------------------------
function [Dtrans,Dreact,defect] = density_energy_diagnostics( ...
    nold,pnew,Mmob,Eold,Enew,par,h)
% Terms in the discrete scalar energy inequality.

if length(nold) > 1
    edge_n = 0.5*(nold(1:end-1)+nold(2:end));
    gradp = diff(pnew)/h;
    Dtrans = par.dt*h*sum(edge_n.*gradp.^2);
else
    Dtrans = 0;
end
Dreact = par.dt*h*sum(Mmob.*(pnew-par.pH).^2);
defect = Enew+Dtrans+Dreact-Eold;
end

%% ------------------------------------------------------------------------
function y = roundoff_project_interval(y,lo,hi,roundoff_tol,label)
low_round = y < lo & y >= lo-roundoff_tol;
high_round = y > hi & y <= hi+roundoff_tol;
y(low_round) = lo;
y(high_round) = hi;

if any(y < lo | y > hi)
    warning('%s exceeds the allowed interval: min=%.3e, max=%.3e.', ...
        label,min(y),max(y));
end
end

%% ------------------------------------------------------------------------
function [n_hist,p_hist,phi1_hist,phi2_hist,n1_hist,n2_hist] = ...
    store_snapshot(sid,n,p,phi1,phi2,n1,n2, ...
    n_hist,p_hist,phi1_hist,phi2_hist,n1_hist,n2_hist)
n_hist(:,sid) = n;
p_hist(:,sid) = p;
phi1_hist(:,sid) = phi1;
phi2_hist(:,sid) = phi2;
n1_hist(:,sid) = n1;
n2_hist(:,sid) = n2;
end

%% ------------------------------------------------------------------------
function [mass,mass1,mass2,fraction1,fraction2] = ...
    compute_masses(n,n1,n2,h)
mass = h*sum(n);
mass1 = h*sum(n1);
mass2 = h*sum(n2);
if mass > 0
    fraction1 = mass1/mass;
    fraction2 = mass2/mass;
else
    fraction1 = NaN;
    fraction2 = NaN;
end
end

%% ------------------------------------------------------------------------
function [core_phi,front_phi] = ...
    compute_core_front_averages(n,p,phi,par)
valid = isfinite(phi);
core_mask = valid & p > par.core_pressure_threshold;
front_mask = valid & p > par.front_pressure_low ...
    & p < par.front_pressure_high;

core_den = sum(n(core_mask));
if core_den > 0
    core_phi = sum(n(core_mask).*phi(core_mask))/core_den;
else
    core_phi = NaN;
end

front_den = sum(n(front_mask));
if front_den > 0
    front_phi = sum(n(front_mask).*phi(front_mask))/front_den;
else
    front_phi = NaN;
end
end

%% ------------------------------------------------------------------------
function E = energy_h(n,par,h)
E = h*sum((1/(par.gamma+1))*n.^(par.gamma+1)-par.pH*n);
end

%% ------------------------------------------------------------------------
function print_case_summary(result,par)
max_energy_jump = max([0;diff(result.energy_hist)]);
energy_decreasing = max_energy_jump <= par.energy_ineq_tol;

[min_phi1,max_phi1] = finite_minmax(result.final_phi1);
[min_phi2,max_phi2] = finite_minmax(result.final_phi2);

fprintf('\nExample 2 summary:\n');
fprintf('  final min/max(n): %.6e / %.6e\n', ...
    min(result.final_n),max(result.final_n));
fprintf('  final min/max(p): %.6e / %.6e\n', ...
    min(result.final_p),max(result.final_p));
fprintf('  final occupied min/max(phi1): %.6e / %.6e\n', ...
    min_phi1,max_phi1);
fprintf('  final occupied min/max(phi2): %.6e / %.6e\n', ...
    min_phi2,max_phi2);
fprintf('  final center phi1/phi2: %.12e / %.12e\n', ...
    result.center_phi1_hist(end),result.center_phi2_hist(end));
fprintf('  final core phi1/phi2: %.12e / %.12e\n', ...
    result.core_phi1_hist(end),result.core_phi2_hist(end));
fprintf('  final front phi1/phi2: %.12e / %.12e\n', ...
    result.front_phi1_hist(end),result.front_phi2_hist(end));
fprintf('  final masses: total=%.12e, mass1=%.12e, mass2=%.12e\n', ...
    result.mass_hist(end),result.mass1_hist(end),result.mass2_hist(end));
fprintf('  final fractions: F1=%.12e, F2=%.12e\n', ...
    result.fraction1_hist(end),result.fraction2_hist(end));
fprintf('  energy monotone: %d, max positive jump=%.3e\n', ...
    energy_decreasing,max_energy_jump);
fprintf('  max discrete energy defect=%.3e\n', ...
    max(result.energy_ineq_defect_hist));
fprintf('  Newton iterations: min=%d, max=%d, mean=%.2f\n', ...
    min(result.newton_it_hist),max(result.newton_it_hist), ...
    mean(result.newton_it_hist));
fprintf('  max Newton residual=%.3e\n',max(result.newton_res_hist));
fprintf('  max scalar source/balance defect=%.3e\n', ...
    max(result.source_balance_defect_hist));
fprintf('  max reaction-allocation error=%.3e\n', ...
    max(result.reaction_allocation_err_hist));
fprintf('  max local reaction overshoot=%.3e\n', ...
    max(result.reaction_fraction_overshoot_hist));
fprintf('  remap overshoot count=%d, max overshoot=%.3e\n', ...
    result.remap_overshoot_count,result.remap_overshoot_max);
fprintf('  max remap simplex error=%.3e\n', ...
    max(result.remap_simplex_err_hist));
fprintf('  max reconstruction error=%.3e\n', ...
    max(result.remap_reconstruction_err_hist));
fprintf('  max species-remap residual=%.3e\n', ...
    max(result.remap_residual_hist));
end

%% ------------------------------------------------------------------------
function [vmin,vmax] = finite_minmax(v)
v = v(isfinite(v));
if isempty(v)
    vmin = NaN;
    vmax = NaN;
else
    vmin = min(v);
    vmax = max(v);
end
end

%% ------------------------------------------------------------------------
