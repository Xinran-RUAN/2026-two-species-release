function result = ms_tumor1d_solver(par,init)
%MS_TUMOR1D_SOLVER
% Reusable 1D solver for a two-species pressure-driven tumor model:
%   1. scalar semi-implicit Onsager density-pressure step;
%   2. normalized frozen-pressure local reaction allocation;
%   3. implicit conservative upwind species remapping driven by the same
%      physical total-density flux.
%
% Required scalar parameter fields:
%   a,b,N,T,dt,gamma,pH
%
% Required reaction function handles:
%   H1(p), H2(p)          source-wise net growth rates;
%   sigma1(p), sigma2(p)  continuous secant coefficients
%                         H_i(p)/(p-pH);
%   k12(p), k21(p)        switching rates 1->2 and 2->1.
%
% Required initial-data function handles:
%   init.n0(x), init.phi1(x)
%
% For backward compatibility, if H1,H2,sigma1,sigma2,k12,k21 are absent
% but scalar fields a1,a2 are supplied, the solver recovers the old
% diagonal affine model H_i(p)=a_i(pH-p), k12=k21=0.
%
% The normalized 2-by-2 matrix exponential is evaluated by a vectorized
% closed formula. No cellwise call to EXPM is used.

par = local_defaults(par);
validate_inputs(par,init);

N = par.N;
h = (par.b-par.a)/N;
x = par.a + ((1:N)'-0.5)*h;
nt = round(par.T/par.dt);
if nt < 1
    error('The requested time step is larger than the final time.');
end
dt = par.T/nt;
par.dt = dt;
time = (0:nt)'*dt;

n0 = init.n0(x);
phi10 = init.phi1(x);
n0 = min(max(n0,0),par.nH);
phi10 = min(max(phi10,0),1);

n1 = n0.*phi10;
n2 = n0-n1;
n = n1+n2;
p = n.^par.gamma;

% Diagnostics. All arrays include the initial time unless otherwise noted.
diag.time = time;
diag.energy = zeros(nt+1,1);
diag.mass = zeros(nt+1,1);
diag.mass1 = zeros(nt+1,1);
diag.mass2 = zeros(nt+1,1);
diag.min_n1 = zeros(nt+1,1);
diag.min_n2 = zeros(nt+1,1);
diag.min_n = zeros(nt+1,1);
diag.max_n = zeros(nt+1,1);
diag.max_p = zeros(nt+1,1);
diag.reconstruction_error = zeros(nt+1,1);
diag.simplex_error = zeros(nt+1,1);

diag.newton_iterations = zeros(nt,1);
diag.newton_residual = zeros(nt,1);
diag.energy_defect = zeros(nt,1);
diag.transport_dissipation = zeros(nt,1);
diag.reaction_dissipation = zeros(nt,1);
diag.scalar_balance_error = zeros(nt,1);
diag.total_mass_balance_error = zeros(nt,1);
diag.allocation_error = zeros(nt,1);
diag.reaction_map_raw_negativity = zeros(nt,1);
diag.remap_residual = zeros(nt,1);
diag.remap_compatibility_error = zeros(nt,1);
diag.remap_raw_negative_rhs = zeros(nt,1);
diag.active_cells = zeros(nt,1);

if par.store_solution_history
    if isempty(par.snapshot_times)
        snapshot_steps = (0:nt)';
    else
        snapshot_steps = unique(round(par.snapshot_times(:)/dt));
        snapshot_steps = max(0,min(nt,snapshot_steps));
        snapshot_steps = unique([0;snapshot_steps;nt]);
    end
    snapshot_times = snapshot_steps*dt;
    par.snapshot_times = snapshot_times;
    nsnapshot = numel(snapshot_steps);

    history.time = snapshot_times;
    history.n = zeros(N,nsnapshot);
    history.n1 = zeros(N,nsnapshot);
    history.n2 = zeros(N,nsnapshot);
    history.p = zeros(N,nsnapshot);
    history.phi1 = NaN(N,nsnapshot);
    history.phi2 = NaN(N,nsnapshot);

    snapshot_id = 1;
    history = store_solution_snapshot( ...
        history,snapshot_id,n1,n2,n,p,par.occupied_tol);
    snapshot_id = snapshot_id+1;
end

[diag.energy(1),diag.mass(1),diag.mass1(1),diag.mass2(1), ...
 diag.min_n1(1),diag.min_n2(1),diag.min_n(1),diag.max_n(1), ...
 diag.max_p(1),diag.reconstruction_error(1),diag.simplex_error(1)] = ...
    state_diagnostics(n1,n2,n,p,par,h);

for k = 1:nt
    n1old = n1;
    n2old = n2;
    nold = n1old+n2old;
    pold = nold.^par.gamma;
    Eold = diag.energy(k);

    % Explicit composition-dependent reaction mobility:
    % M^k = -sum_i n_i^k sigma_i(p^k).
    sigma1old = par.sigma1(pold);
    sigma2old = par.sigma2(pold);
    Mmob = -(n1old.*sigma1old+n2old.*sigma2old);

    mobility_scale = max([1;abs(Mmob(:))]);
    mobility_tol = max(par.roundoff_tol*mobility_scale, ...
        100*eps(mobility_scale)*mobility_scale);
    if min(Mmob) < -mobility_tol
        error('Negative reaction mobility: min(M)=%.3e.',min(Mmob));
    end
    Mmob(Mmob<0) = 0;

    [nnew,pnew,J,newton_info] = density_step(nold,Mmob,par,h);
    Snew = Mmob.*(par.pH-pnew);
    nhat = nold+dt*Snew;

    scalar_balance = nnew+(dt/h)*(J(2:end)-J(1:end-1))-nhat;
    scalar_balance_error = max(abs(scalar_balance));

    [nhat_species,allocation_info] = local_reaction_allocation( ...
        n1old,n2old,nhat,pnew,par);

    [n_species,theta,remap_info] = conservative_species_remap( ...
        nnew,nhat_species,J,pnew,dt,h,par);

    n1 = n_species(:,1);
    n2 = n_species(:,2);
    n = n1+n2;
    p = n.^par.gamma;

    Enew = energy_h(n,par,h);
    [Dtrans,Dreact,Edefect] = energy_diagnostics( ...
        nold,pnew,Mmob,Eold,Enew,par,h);

    diag.newton_iterations(k) = newton_info.iter;
    diag.newton_residual(k) = newton_info.residual;
    diag.energy_defect(k) = Edefect;
    diag.transport_dissipation(k) = Dtrans;
    diag.reaction_dissipation(k) = Dreact;
    diag.scalar_balance_error(k) = scalar_balance_error;
    diag.total_mass_balance_error(k) = abs( ...
        h*sum(nnew-nold)-dt*h*sum(Snew));
    diag.allocation_error(k) = allocation_info.allocation_error;
    diag.reaction_map_raw_negativity(k) = ...
        allocation_info.max_raw_negative_component;
    diag.remap_residual(k) = remap_info.species_residual;
    diag.remap_compatibility_error(k) = ...
        remap_info.mass_compatibility_defect;
    diag.remap_raw_negative_rhs(k) = remap_info.max_raw_negative_rhs;
    diag.active_cells(k) = remap_info.active_cell_count;

    [diag.energy(k+1),diag.mass(k+1),diag.mass1(k+1),diag.mass2(k+1), ...
     diag.min_n1(k+1),diag.min_n2(k+1),diag.min_n(k+1),diag.max_n(k+1), ...
     diag.max_p(k+1),diag.reconstruction_error(k+1),diag.simplex_error(k+1)] = ...
        state_diagnostics(n1,n2,n,p,par,h,theta);

    if par.store_solution_history && snapshot_id <= nsnapshot && ...
            k == snapshot_steps(snapshot_id)
        history = store_solution_snapshot( ...
            history,snapshot_id,n1,n2,n,p,par.occupied_tol);
        snapshot_id = snapshot_id+1;
    end

    if par.progress_every > 0 && ...
            (mod(k,par.progress_every)==0 || k==nt)
        fprintf(['step %6d/%6d, t=%8.5f, Newton=%2d, res=%.3e, ', ...
            'Edef=%.3e, remap=%.3e, recon=%.3e\n'], ...
            k,nt,k*dt,newton_info.iter,newton_info.residual,Edefect, ...
            remap_info.species_residual,diag.reconstruction_error(k+1));
    end
end

if par.store_solution_history && snapshot_id ~= nsnapshot+1
    error('The solver did not store every requested solution snapshot.');
end

result.par = par;
result.init = init;
result.x = x;
result.h = h;
result.dt = dt;
result.nt = nt;
result.n1 = n1;
result.n2 = n2;
result.n = n;
result.p = p;
result.phi1 = species_fraction(n1,n2,par.occupied_tol);
result.phi2 = species_fraction(n2,n1,par.occupied_tol);
result.diagnostics = diag;
if par.store_solution_history
    result.history = history;
end
end

%% ========================================================================
function history = store_solution_snapshot( ...
    history,snapshot_id,n1,n2,n,p,occupied_tol)
history.n(:,snapshot_id) = n;
history.n1(:,snapshot_id) = n1;
history.n2(:,snapshot_id) = n2;
history.p(:,snapshot_id) = p;
history.phi1(:,snapshot_id) = species_fraction(n1,n2,occupied_tol);
history.phi2(:,snapshot_id) = species_fraction(n2,n1,occupied_tol);
end

%% ========================================================================
function par = local_defaults(par)
% Backward-compatible conversion of the former diagonal affine model.
if ~isfield(par,'H1') && isfield(par,'a1')
    a1 = par.a1;
    pH = par.pH;
    par.H1 = @(p) a1.*(pH-p);
end
if ~isfield(par,'H2') && isfield(par,'a2')
    a2 = par.a2;
    pH = par.pH;
    par.H2 = @(p) a2.*(pH-p);
end
if ~isfield(par,'sigma1') && isfield(par,'a1')
    a1 = par.a1;
    par.sigma1 = @(p) -a1+zeros(size(p));
end
if ~isfield(par,'sigma2') && isfield(par,'a2')
    a2 = par.a2;
    par.sigma2 = @(p) -a2+zeros(size(p));
end
if ~isfield(par,'k12')
    par.k12 = @(p) zeros(size(p));
end
if ~isfield(par,'k21')
    par.k21 = @(p) zeros(size(p));
end

if ~isfield(par,'nH'), par.nH = par.pH^(1/par.gamma); end
if ~isfield(par,'newton_tol'), par.newton_tol = 1e-11; end
if ~isfield(par,'newton_maxit'), par.newton_maxit = 50; end
if ~isfield(par,'line_search_maxit'), par.line_search_maxit = 30; end
if ~isfield(par,'armijo'), par.armijo = 1e-4; end
if ~isfield(par,'warn_residual_tol'), par.warn_residual_tol = 1e-8; end
if ~isfield(par,'abort_on_newton_failure')
    par.abort_on_newton_failure = true;
end
if ~isfield(par,'roundoff_tol'), par.roundoff_tol = 1e-12; end
if ~isfield(par,'remap_mass_tol'), par.remap_mass_tol = 1e-13; end
if ~isfield(par,'remap_negative_tol'), par.remap_negative_tol = 1e-11; end
if ~isfield(par,'occupied_tol'), par.occupied_tol = 1e-12; end
if ~isfield(par,'progress_every'), par.progress_every = 0; end
if ~isfield(par,'store_solution_history')
    par.store_solution_history = false;
end
if ~isfield(par,'snapshot_times'), par.snapshot_times = []; end
end

function validate_inputs(par,init)
required_scalar = {'a','b','N','T','dt','gamma','pH'};
for q = 1:numel(required_scalar)
    if ~isfield(par,required_scalar{q})
        error('Missing parameter par.%s.',required_scalar{q});
    end
end

required_reaction = {'H1','H2','sigma1','sigma2','k12','k21'};
for q = 1:numel(required_reaction)
    name = required_reaction{q};
    if ~isfield(par,name) || ~isa(par.(name),'function_handle')
        error('par.%s must be a function handle.',name);
    end
end

if ~isfield(init,'n0') || ~isa(init.n0,'function_handle')
    error('init.n0 must be a function handle.');
end
if ~isfield(init,'phi1') || ~isa(init.phi1,'function_handle')
    error('init.phi1 must be a function handle.');
end
if par.N < 2 || par.N~=round(par.N) || ...
        par.b <= par.a || par.T <= 0 || par.dt <= 0
    error('Invalid grid or time parameters.');
end
if par.gamma <= 1 || par.pH <= 0
    error('Require gamma>1 and pH>0.');
end

expected_nH = par.pH^(1/par.gamma);
if abs(par.nH-expected_nH) > 1e-12*(1+expected_nH)
    error('par.nH must equal par.pH^(1/par.gamma).');
end

% Check the reaction data on the invariant pressure interval.
ptest = linspace(0,par.pH,65)';
H1 = evaluate_reaction_handle(par.H1,ptest,'H1');
H2 = evaluate_reaction_handle(par.H2,ptest,'H2');
sigma1 = evaluate_reaction_handle(par.sigma1,ptest,'sigma1');
sigma2 = evaluate_reaction_handle(par.sigma2,ptest,'sigma2');
k12 = evaluate_reaction_handle(par.k12,ptest,'k12');
k21 = evaluate_reaction_handle(par.k21,ptest,'k21');

scale = max([1;abs(H1);abs(H2);abs(sigma1);abs(sigma2);abs(k12);abs(k21)]);
tol = 1e-10*scale;

if abs(H1(end)) > tol || abs(H2(end)) > tol
    error('H1 and H2 must vanish at p=pH.');
end
if any(H1 < -tol) || any(H2 < -tol)
    error('H1 and H2 must be nonnegative on [0,pH].');
end
if any(sigma1 > tol) || any(sigma2 > tol)
    error('sigma1 and sigma2 must be nonpositive on [0,pH].');
end
if any(k12 < -tol) || any(k21 < -tol)
    error('k12 and k21 must be nonnegative on [0,pH].');
end

secant_defect = max([ ...
    max(abs(H1-sigma1.*(ptest-par.pH))), ...
    max(abs(H2-sigma2.*(ptest-par.pH)))]);
if secant_defect > 1e-9*scale
    error('The supplied H_i and sigma_i are inconsistent.');
end
end

function value = evaluate_reaction_handle(fun,p,name)
value = fun(p);
if ~isnumeric(value) || ~isequal(size(value),size(p)) || ...
        any(~isfinite(value))
    error('par.%s must return a finite array with the same size as p.',name);
end
value = value(:);
end

%% ========================================================================
function [nnew,pnew,J,info] = density_step(nold,Mmob,par,h)
N = length(nold);
dt = par.dt;
n = min(max(nold,0),par.nH);

edge_n = zeros(N+1,1);
edge_n(2:N) = 0.5*(nold(1:N-1)+nold(2:N));

info.converged = false;
info.iter = 0;
info.residual = Inf;
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
        error('Density Newton produced a non-finite direction.');
    end

    lambda = 1;
    accepted = false;
    res2 = norm(R,2);
    best_n = n;
    best_res = res2;

    for ls = 1:par.line_search_maxit
        ntrial = min(max(n+lambda*dn,0),par.nH);
        Rtrial = density_residual(ntrial,nold,edge_n,Mmob,par,h);
        trial_res = norm(Rtrial,2);
        if isfinite(trial_res) && trial_res < best_res
            best_res = trial_res;
            best_n = ntrial;
        end
        if isfinite(trial_res) && ...
                trial_res <= (1-par.armijo*lambda)*res2
            accepted = true;
            break;
        end
        lambda = 0.5*lambda;
    end

    if accepted
        n = ntrial;
    elseif best_res < res2
        n = best_n;
    else
        info.iter = it;
        info.residual = res_inf;
        break;
    end
    info.iter = it;
end

nnew = min(max(n,0),par.nH);
pnew = nnew.^par.gamma;
Rfinal = density_residual(nnew,nold,edge_n,Mmob,par,h);
info.residual = norm(Rfinal,inf);
if info.residual <= stop_tol
    info.converged = true;
end
if ~info.converged && info.residual > par.warn_residual_tol
    msg = sprintf( ...
        'Density Newton did not converge: residual=%.3e.',info.residual);
    if par.abort_on_newton_failure
        error(msg);
    else
        warning('%s',msg);
    end
end

F = zeros(N+1,1);
F(2:N) = edge_n(2:N).*(pnew(2:N)-pnew(1:N-1))/h;
J = -F;
end

function R = density_residual(n,nold,edge_n,Mmob,par,h)
N = length(n);
p = n.^par.gamma;
F = zeros(N+1,1);
F(2:N) = edge_n(2:N).*(p(2:N)-p(1:N-1))/h;
dhF = (F(2:N+1)-F(1:N))/h;
R = n-nold-par.dt*dhF+par.dt*Mmob.*(p-par.pH);
end

function [R,Jmat] = density_residual_jac( ...
    n,nold,edge_n,Mmob,par,h)
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
    diag_entry = 1+par.dt*Mmob(j)*dpdn(j) ...
        +coef*(Aleft+Aright)*dpdn(j);

    cnt = cnt+1;
    row(cnt) = j;
    col(cnt) = j;
    val(cnt) = diag_entry;

    if j>1
        cnt = cnt+1;
        row(cnt) = j;
        col(cnt) = j-1;
        val(cnt) = -coef*Aleft*dpdn(j-1);
    end
    if j<N
        cnt = cnt+1;
        row(cnt) = j;
        col(cnt) = j+1;
        val(cnt) = -coef*Aright*dpdn(j+1);
    end
end
Jmat = sparse(row(1:cnt),col(1:cnt),val(1:cnt),N,N);
end

%% ========================================================================
function [nhat_species,info] = local_reaction_allocation( ...
    n1old,n2old,nhat,pnew,par)
% Normalized frozen-pressure map
%
%   Psi_dt(phi,p) = exp(dt*A(p))*phi
%                   --------------------------------,
%                   1^T exp(dt*A(p))*phi
%
% where A=R^T. For
%
%   R = [H1-k12, k12;
%        k21,    H2-k21],
%
% the matrix used here is
%
%   A = [H1-k12, k21;
%        k12,    H2-k21].
%
% A vectorized closed formula is used. Writing
%
%   A = m*I + B,  B=[q,b;c,-q],  delta^2=q^2+b*c,
%
% the common positive factor exp(dt*m)*cosh(dt*delta) cancels in the
% normalization. Hence only
%
%   I + tanh(dt*delta)/delta * B
%
% needs to be applied. The coefficient has the limit dt when delta=0.

old_total = n1old+n2old;
active = old_total > par.occupied_tol;
frac1 = zeros(size(old_total));
max_raw_negative = 0;

nhat_scale = max([1;abs(nhat(:))]);
nhat_tol = max(par.roundoff_tol*nhat_scale, ...
    100*eps(nhat_scale)*nhat_scale);
if min(nhat) < -nhat_tol
    error('Negative source-updated total density: min(nhat)=%.3e.', ...
        min(nhat));
end
nhat(nhat<0) = 0;

if any(active)
    p = pnew(active);
    H1 = par.H1(p);
    H2 = par.H2(p);
    k12 = par.k12(p);
    k21 = par.k21(p);

    reaction_scale = max([1;abs(H1);abs(H2);abs(k12);abs(k21)]);
    reaction_tol = max(par.roundoff_tol*reaction_scale, ...
        100*eps(reaction_scale)*reaction_scale);
    if min(k12) < -reaction_tol || min(k21) < -reaction_tol
        error('A switching rate became negative in the reaction map.');
    end
    k12(k12<0) = 0;
    k21(k21<0) = 0;

    % Entries of A=R^T.
    A11 = H1-k12;
    A12 = k21;
    A21 = k12;
    A22 = H2-k21;

    q = 0.5*(A11-A22);
    delta2 = q.^2+A12.*A21;
    delta = sqrt(max(delta2,0));

    coeff = par.dt*ones(size(delta));
    nonsingular = delta > sqrt(eps);
    coeff(nonsingular) = ...
        tanh(par.dt*delta(nonsingular))./delta(nonsingular);

    v1 = n1old(active);
    v2 = n2old(active);

    % This vector is proportional to exp(dt*A)*[v1;v2].
    y1 = v1+coeff.*(q.*v1+A12.*v2);
    y2 = v2+coeff.*(A21.*v1-q.*v2);

    map_scale = max([1;abs(v1);abs(v2);abs(y1);abs(y2)]);
    map_tol = max(par.roundoff_tol*map_scale, ...
        100*eps(map_scale)*map_scale);
    max_raw_negative = max(0,-min([y1;y2]));
    if min([y1;y2]) < -map_tol
        error(['The normalized reaction map produced a negative ', ...
            'component: min(y)=%.3e.'],min([y1;y2]));
    end
    y1(y1<0) = 0;
    y2(y2<0) = 0;

    den = y1+y2;
    if any(~isfinite(den)) || any(den<=0)
        error('The normalized reaction map has a nonpositive denominator.');
    end

    tmp = y1./den;
    frac1(active) = min(max(tmp,0),1);
end

nhat1 = nhat.*frac1;
nhat2 = nhat-nhat1;
nhat_species = [nhat1,nhat2];

info.allocation_error = max(abs(sum(nhat_species,2)-nhat));
info.max_raw_negative_component = max_raw_negative;
end

%% ========================================================================
function [n_species,theta,info] = conservative_species_remap( ...
    nnew,nhat_species,J,pnew,dt,h,par)
% Pressure-ordered solution of the implicit conservative upwind remapping.

N = length(nnew);
lambda = dt/h;
Jplus = max(J,0);
Jminus = min(J,0);

d = nnew+lambda*(Jplus(2:N+1)-Jminus(1:N));
mass_scale = max([1;abs(nnew(:));abs(nhat_species(:));lambda*abs(J(:))]);
mass_tol = max(par.remap_mass_tol*mass_scale, ...
    100*eps(mass_scale)*mass_scale);
neg_tol = max(par.remap_negative_tol*mass_scale, ...
    100*eps(mass_scale)*mass_scale);

if min(d) < -neg_tol
    error('Negative remapping diagonal: min(d)=%.3e.',min(d));
end
d(d<0) = 0;
active = d > mass_tol;

[~,order] = sort(pnew,'descend');
theta = zeros(N,2);
processed = false(N,1);
max_compat = 0;
max_raw_negative = 0;

for q = 1:N
    j = order(q);
    if ~active(j), continue; end

    rhs = nhat_species(j,:);

    % Left inflow through edge j, between cells j-1 and j.
    if j>1 && J(j)>0
        if ~processed(j-1)
            error('Pressure ordering failed at left inflow into cell %d.',j);
        end
        rhs = rhs+lambda*J(j)*theta(j-1,:);
    end

    % Right inflow through edge j+1, between cells j and j+1.
    if j<N && J(j+1)<0
        if ~processed(j+1)
            error('Pressure ordering failed at right inflow into cell %d.',j);
        end
        rhs = rhs-lambda*J(j+1)*theta(j+1,:);
    end

    max_raw_negative = max(max_raw_negative,max(0,-min(rhs)));
    if min(rhs) < -neg_tol
        error('Negative species remapping RHS at cell %d: %.3e.', ...
            j,min(rhs));
    end
    rhs(rhs<0) = 0;
    rhs_mass = sum(rhs);
    max_compat = max(max_compat,abs(rhs_mass-d(j)));

    if rhs_mass > mass_tol
        theta(j,:) = rhs/rhs_mass;
    else
        local_mass = sum(nhat_species(j,:));
        if local_mass > mass_tol
            theta(j,:) = nhat_species(j,:)/local_mass;
        else
            theta(j,:) = [0.5,0.5];
        end
    end

    theta(j,:) = max(theta(j,:),0);
    theta(j,:) = theta(j,:)/sum(theta(j,:));
    processed(j) = true;
end

% Numerically inactive cells carry no resolved outgoing mass. Assign a
% harmless local fraction only when a tiny positive nnew must be rebuilt.
for j = 1:N
    if active(j), continue; end
    if nnew(j)>0
        local_mass = sum(nhat_species(j,:));
        if local_mass>mass_tol
            theta(j,:) = nhat_species(j,:)/local_mass;
        else
            theta(j,:) = [0.5,0.5];
        end
    end
end

n1new = nnew.*theta(:,1);
n2new = nnew-n1new;
n_species = [n1new,n2new];

Jspecies = zeros(N+1,2);
for e = 2:N
    if J(e)>0
        Jspecies(e,:) = J(e)*theta(e-1,:);
    elseif J(e)<0
        Jspecies(e,:) = J(e)*theta(e,:);
    end
end

Rspecies = n_species ...
    +lambda*(Jspecies(2:N+1,:)-Jspecies(1:N,:)) ...
    -nhat_species;

resolved = active | (nnew>par.occupied_tol);
if any(resolved)
    simplex_error = max(abs(sum(theta(resolved,:),2)-1));
else
    simplex_error = 0;
end

info.active_cell_count = nnz(active);
info.simplex_error = simplex_error;
info.reconstruction_error = max(abs(sum(n_species,2)-nnew));
info.species_residual = max(abs(Rspecies(:)));
info.mass_compatibility_defect = max_compat;
info.max_raw_negative_rhs = max_raw_negative;
end

%% ========================================================================
function [E,mass,mass1,mass2,minn1,minn2,minn,maxn,maxp,recon,simplex] = ...
    state_diagnostics(n1,n2,n,p,par,h,theta)
if nargin<7
    theta = [];
end

E = energy_h(n,par,h);
mass = h*sum(n);
mass1 = h*sum(n1);
mass2 = h*sum(n2);
minn1 = min(n1);
minn2 = min(n2);
minn = min(n);
maxn = max(n);
maxp = max(p);
recon = max(abs(n1+n2-n));

if isempty(theta)
    occupied = n>par.occupied_tol;
    if any(occupied)
        phi1 = n1(occupied)./n(occupied);
        phi2 = n2(occupied)./n(occupied);
        simplex = max(abs(phi1+phi2-1));
    else
        simplex = 0;
    end
else
    resolved = n>par.occupied_tol;
    if any(resolved)
        simplex = max(abs(sum(theta(resolved,:),2)-1));
    else
        simplex = 0;
    end
end
end

function [Dtrans,Dreact,defect] = energy_diagnostics( ...
    nold,pnew,Mmob,Eold,Enew,par,h)
edge_n = 0.5*(nold(1:end-1)+nold(2:end));
gradp = diff(pnew)/h;
Dtrans = par.dt*h*sum(edge_n.*gradp.^2);
Dreact = par.dt*h*sum(Mmob.*(pnew-par.pH).^2);
defect = Enew+Dtrans+Dreact-Eold;
end

function E = energy_h(n,par,h)
E = h*sum(n.^(par.gamma+1)/(par.gamma+1)-par.pH*n);
end

function phi1 = species_fraction(n1,n2,tol)
n = n1+n2;
phi1 = NaN(size(n));
idx = n>tol;
phi1(idx) = n1(idx)./n(idx);
end
