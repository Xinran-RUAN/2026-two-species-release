function run_section4_4_stiff_pressure(outputDir)
%RUN_SECTION4_4_STIFF_PRESSURE
% Section 4.4: fixed-grid stiff-pressure limit and moving species interface.
%
% This program performs the numerical runs and saves the final fields and
% diagnostics.  It uses the current three-step scheme:
%
%   Step 1. Semi-implicit scalar density-pressure update with frozen
%           mobility M^k=n^k.
%   Step 2. Normalized frozen-pressure reaction allocation.
%           Here R(p)=(1-p)I, so the normalized reaction map leaves the
%           local composition unchanged.
%   Step 3. Implicit conservative species remapping driven by exactly the
%           same total physical mass flux as the scalar step.
%
% The separate files
%   plot_section4_4_density_pressure.m
%   make_section4_4_tables.m
% produce the paper figure and the two LaTeX tables.
%
% The full production run uses N=6400, dt=1e-4, T=1 and seven gamma values.
% It may take a substantial amount of time.  A checkpoint is written after
% every completed gamma value, so an interrupted run can be resumed.

clc; close all;
if nargin < 1 || isempty(outputDir), outputDir=fullfile(tempdir,'two_species_rerun','stiff_pressure'); end

%% ========================================================================
% Parameters
% ========================================================================

par=model.stiff_pressure_parameters;

fprintf('Code version: %s\n',par.code_version);

%% ========================================================================
% Grid, time, output
% ========================================================================

N = par.N;
h = (par.b-par.a)/N;
x = par.a + ((1:N)'-0.5)*h;

nt = round(par.T/par.dt);
par.dt = par.T/nt;
time_grid = (0:nt)'*par.dt;


% outputDir is supplied by the caller; numerical parameters are unchanged.
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

dataFile = fullfile(outputDir,'section4_4_stiff_pressure_data.mat');

%% ========================================================================
% Exact Hele--Shaw reference
% ========================================================================

reference.R_t = exact_outer_boundary(time_grid,par.R0);
reference.X_t = exact_internal_interface( ...
    reference.R_t,par.R0,par.X0);

reference.R_T = reference.R_t(end);
reference.X_T = reference.X_t(end);

reference.n_T = double(abs(x) <= reference.R_T);
reference.p_T = max( ...
    1-cosh(x)./cosh(reference.R_T), ...
    0);

%% ========================================================================
% Allocate or resume
% ========================================================================

gammas = par.gammas(:).';
ng = numel(gammas);

final_n = NaN(N,ng);
final_p = NaN(N,ng);
final_n1 = NaN(N,ng);
final_n2 = NaN(N,ng);

error_n = NaN(1,ng);
error_p = NaN(1,ng);
error_R = NaN(1,ng);
error_X = NaN(1,ng);

R_num = NaN(1,ng);
X_num = NaN(1,ng);

complementarity_defect = NaN(1,ng);
interior_pressure_residual = NaN(1,ng);

newton_mean_iterations = NaN(1,ng);
newton_max_iterations = NaN(1,ng);
max_newton_residual = NaN(1,ng);
max_reconstruction_error = NaN(1,ng);
max_simplex_error = NaN(1,ng);
max_remap_residual = NaN(1,ng);

completed = false(1,ng);

if par.resume && exist(dataFile,'file') == 2
    old = load(dataFile);

    compatible = isfield(old,'par') ...
        && isfield(old,'x') ...
        && isequal(old.par.gammas(:).',gammas) ...
        && old.par.N == par.N ...
        && abs(old.par.dt-par.dt) <= 10*eps ...
        && abs(old.par.T-par.T) <= 10*eps ...
        && numel(old.x) == numel(x) ...
        && max(abs(old.x(:)-x(:))) <= 10*eps;

    if compatible
        fprintf('Resuming from %s\n',dataFile);

        names = { ...
            'final_n','final_p','final_n1','final_n2', ...
            'error_n','error_p','error_R','error_X', ...
            'R_num','X_num', ...
            'complementarity_defect','interior_pressure_residual', ...
            'newton_mean_iterations','newton_max_iterations', ...
            'max_newton_residual','max_reconstruction_error', ...
            'max_simplex_error','max_remap_residual','completed'};

        for q = 1:numel(names)
            name = names{q};
            if isfield(old,name)
                eval([name ' = old.(name);']); %#ok<EVLDIR>
            end
        end
    else
        warning('Existing checkpoint is incompatible and will be overwritten.');
    end
end

%% ========================================================================
% Production runs
% ========================================================================

for ig = 1:ng
    gamma = gammas(ig);

    if completed(ig)
        fprintf('\nSkipping gamma=%g: checkpoint already completed.\n',gamma);
        continue;
    end

    fprintf('\n============================================================\n');
    fprintf('Running gamma=%g, N=%d, dt=%.3e, T=%.3g\n', ...
        gamma,N,par.dt,par.T);
    fprintf('============================================================\n');

    % Initial pressure and gamma-dependent initial total density.
    [p0,n0,n1_0,n2_0]=model.stiff_initial(x,gamma,par);

    result = run_one_gamma( ...
        n1_0,n2_0,x,h,nt,gamma,par);

    final_n(:,ig) = result.n;
    final_p(:,ig) = result.p;
    final_n1(:,ig) = result.n1;
    final_n2(:,ig) = result.n2;

    error_n(ig) = discrete_l1( ...
        result.n-reference.n_T,h);

    error_p(ig) = discrete_l1( ...
        result.p-reference.p_T,h);

    R_num(ig) = 0.5*h*sum(result.n);
    X_num(ig) = 0.5*h*sum(result.n1-result.n2);

    error_R(ig) = abs(R_num(ig)-reference.R_T);
    error_X(ig) = abs(X_num(ig)-reference.X_T);

    complementarity_defect(ig) = ...
        h*sum(abs(result.p.*(1-result.n)));

    interior_pressure_residual(ig) = ...
        compute_interior_pressure_residual( ...
        result.p,reference.p_T,h, ...
        par.interior_pressure_cutoff);

    newton_mean_iterations(ig) = result.newton_mean_iterations;
    newton_max_iterations(ig) = result.newton_max_iterations;
    max_newton_residual(ig) = result.max_newton_residual;
    max_reconstruction_error(ig) = result.max_reconstruction_error;
    max_simplex_error(ig) = result.max_simplex_error;
    max_remap_residual(ig) = result.max_remap_residual;

    completed(ig) = true;

    fprintf('\nFinal diagnostics for gamma=%g\n',gamma);
    fprintf('  e_n        = %.12e\n',error_n(ig));
    fprintf('  e_p        = %.12e\n',error_p(ig));
    fprintf('  e_R        = %.12e\n',error_R(ig));
    fprintf('  e_X        = %.12e\n',error_X(ig));
    fprintf('  C_gamma    = %.12e\n',complementarity_defect(ig));
    fprintf('  H_gamma    = %.12e\n',interior_pressure_residual(ig));
    fprintf('  R_num      = %.12e\n',R_num(ig));
    fprintf('  X_num      = %.12e\n',X_num(ig));
    fprintf('  Newton mean/max = %.2f / %d\n', ...
        newton_mean_iterations(ig),newton_max_iterations(ig));
    fprintf('  max Newton residual = %.3e\n',max_newton_residual(ig));
    fprintf('  max reconstruction  = %.3e\n',max_reconstruction_error(ig));
    fprintf('  max simplex error   = %.3e\n',max_simplex_error(ig));
    fprintf('  max remap residual  = %.3e\n',max_remap_residual(ig));

    if par.save_after_each_gamma
        save_checkpoint();
    end
end

save_checkpoint();

fprintf('\nSection 4.4 run completed.\n');
fprintf('Saved data:\n  %s\n',dataFile);

%% ========================================================================
% Nested checkpoint writer
% ========================================================================

    function save_checkpoint()
        save(dataFile, ...
            'x','h','time_grid','gammas','par','reference', ...
            'final_n','final_p','final_n1','final_n2', ...
            'error_n','error_p','error_R','error_X', ...
            'R_num','X_num', ...
            'complementarity_defect','interior_pressure_residual', ...
            'newton_mean_iterations','newton_max_iterations', ...
            'max_newton_residual','max_reconstruction_error', ...
            'max_simplex_error','max_remap_residual', ...
            'completed','-v7.3');
    end

end

%% ========================================================================
% One gamma run
% ========================================================================

function result = run_one_gamma(n1_0,n2_0,x,h,nt,gamma,par)

N = numel(x);

n1 = max(n1_0,0);
n2 = max(n2_0,0);
n = n1+n2;
p = n.^gamma;

newtonIterations = zeros(nt,1);
maxNewtonResidual = 0;
maxReconstructionError = 0;
maxSimplexError = 0;
maxRemapResidual = 0;

for k = 1:nt
    nold = n;
    n1old = n1;
    n2old = n2;

    oldReconstruction = max(abs(n1old+n2old-nold));
    if oldReconstruction > par.reconstruction_tol
        error('Reconstruction failed before step %d: %.3e.', ...
            k,oldReconstruction);
    end

    % Step 1.  Since H_1=H_2=1-p, sigma_1=sigma_2=-1 and M^k=n^k.
    Mmob = nold;

    [nnew,pnew,F,info] = solve_density_step( ...
        nold,Mmob,gamma,par,h);

    newtonIterations(k) = info.iter;
    maxNewtonResidual = max(maxNewtonResidual,info.residual);

    % Paper flux F=n_edge*delta_h p has the opposite sign to the physical
    % mass flux J=-n grad p.
    J = -F;

    % Source-updated total density.
    source = Mmob.*(par.pH-pnew);
    sourceTol = par.roundoff_tol*(1+max(abs(source)));

    if any(source < -sourceTol)
        error('Negative source at step %d: %.3e.',k,min(source));
    end
    source(source < 0) = 0;

    nhat = nold+par.dt*source;

    % Step 2.  R(p)=(1-p)I, so exp(dt A) is a scalar multiple of I.
    % The normalized reaction map therefore leaves phi unchanged.
    oldActive = nold > 0;
    phi1old = zeros(N,1);
    phi1old(oldActive) = n1old(oldActive)./nold(oldActive);

    nhat1 = nhat.*phi1old;
    nhat2 = nhat-nhat1;
    nhatSpecies = [nhat1,nhat2];

    % Step 3.  Implicit conservative species remapping.
    % The remapping routine constructs a topological order directly from
    % the actual directed flux graph.  This avoids ambiguous tie-breaking
    % at equal-pressure cells, especially the two cells adjacent to x=0.
    [species,remapInfo] = conservative_species_remap( ...
        nnew,nhatSpecies,J,par.dt,h,par);

    n1 = species(:,1);
    n2 = species(:,2);
    n = nnew;
    p = pnew;

    maxReconstructionError = max( ...
        maxReconstructionError,remapInfo.reconstruction_error);

    maxSimplexError = max( ...
        maxSimplexError,remapInfo.simplex_error);

    maxRemapResidual = max( ...
        maxRemapResidual,remapInfo.species_residual);

    if mod(k,par.progress_every) == 0 || k == nt
        fprintf(['  step %5d/%5d, t=%.4f, Newton=%d, res=%.3e, ', ...
            'remap=%.3e\n'], ...
            k,nt,k*par.dt,info.iter,info.residual, ...
            remapInfo.species_residual);
    end
end

result.n = n;
result.p = p;
result.n1 = n1;
result.n2 = n2;

result.newton_mean_iterations = mean(newtonIterations);
result.newton_max_iterations = max(newtonIterations);
result.max_newton_residual = maxNewtonResidual;
result.max_reconstruction_error = maxReconstructionError;
result.max_simplex_error = maxSimplexError;
result.max_remap_residual = maxRemapResidual;
end

%% ========================================================================
% Scalar density-pressure step
% ========================================================================

function [nnew,pnew,F,info] = solve_density_step( ...
    nold,Mmob,gamma,par,h)

N = numel(nold);

% The invariant interval is [0,n_H]=[0,1] because p_H=1.
n = min(max(nold,0),1);

edgeN = zeros(N+1,1);
if N > 1
    edgeN(2:N) = 0.5*(nold(1:N-1)+nold(2:N));
end

stopTol = par.newton_tol*(1+norm(nold,inf));

info.converged = false;
info.iter = 0;
info.residual = Inf;

for it = 1:par.newton_maxit
    [R,Jmat] = density_residual_jacobian( ...
        n,nold,edgeN,Mmob,gamma,par,h);

    residualInf = norm(R,inf);

    if residualInf <= stopTol
        info.converged = true;
        info.iter = it-1;
        info.residual = residualInf;
        break;
    end

    dn = -(Jmat\R);

    if any(~isfinite(dn))
        error('Non-finite Newton direction.');
    end

    lambda = 1;
    currentNorm = norm(R,2);
    accepted = false;

    bestN = n;
    bestNorm = currentNorm;

    for ls = 1:par.line_search_maxit
        ntrial = min(max(n+lambda*dn,0),1);

        Rtrial = density_residual( ...
            ntrial,nold,edgeN,Mmob,gamma,par,h);

        trialNorm = norm(Rtrial,2);

        if isfinite(trialNorm) && trialNorm < bestNorm
            bestNorm = trialNorm;
            bestN = ntrial;
        end

        if isfinite(trialNorm) ...
                && trialNorm <= (1-par.armijo*lambda)*currentNorm
            accepted = true;
            break;
        end

        lambda = 0.5*lambda;
    end

    if accepted
        n = ntrial;
    elseif bestNorm < currentNorm
        n = bestN;
    else
        warning('Newton line search failed at iteration %d.',it);
        break;
    end

    info.iter = it;
end

nnew = roundoff_project_interval( ...
    n,0,1,par.roundoff_tol,'density');

pnew = nnew.^gamma;

Rfinal = density_residual( ...
    nnew,nold,edgeN,Mmob,gamma,par,h);

info.residual = norm(Rfinal,inf);

if info.residual <= stopTol
    info.converged = true;
end

if ~info.converged && info.residual > par.warn_residual_tol
    message = sprintf( ...
        'Density Newton failed: residual %.3e.',info.residual);

    if par.abort_on_newton_failure
        error(message);
    else
        warning('%s',message);
    end
end

F = zeros(N+1,1);
if N > 1
    F(2:N) = edgeN(2:N).*(pnew(2:N)-pnew(1:N-1))/h;
end
end

function R = density_residual( ...
    n,nold,edgeN,Mmob,gamma,par,h)

N = numel(n);
p = n.^gamma;

F = zeros(N+1,1);
if N > 1
    F(2:N) = edgeN(2:N).*(p(2:N)-p(1:N-1))/h;
end

divF = (F(2:N+1)-F(1:N))/h;

R = n-nold-par.dt*divF ...
    +par.dt*Mmob.*(p-par.pH);
end

function [R,Jmat] = density_residual_jacobian( ...
    n,nold,edgeN,Mmob,gamma,par,h)

N = numel(n);

R = density_residual( ...
    n,nold,edgeN,Mmob,gamma,par,h);

dpdn = gamma*n.^(gamma-1);
coef = par.dt/(h*h);

Aleft = edgeN(1:N);
Aright = edgeN(2:N+1);

mainDiag = 1 ...
    +par.dt*Mmob.*dpdn ...
    +coef*(Aleft+Aright).*dpdn;

lowerDiag = -coef*Aleft(2:N).*dpdn(1:N-1);
upperDiag = -coef*Aright(1:N-1).*dpdn(2:N);

rows = [(1:N)'; (2:N)'; (1:N-1)'];
cols = [(1:N)'; (1:N-1)'; (2:N)'];
values = [mainDiag; lowerDiag; upperDiag];

Jmat = sparse(rows,cols,values,N,N);
end

%% ========================================================================
% Conservative species remapping
% ========================================================================

function [species,info] = conservative_species_remap( ...
    nnew,nhatSpecies,J,dt,h,par)
%CONSERVATIVE_SPECIES_REMAP
% Implicit conservative remapping solved in a topological order determined
% by the actual directed total-mass flux.
%
% In exact arithmetic, every nonzero flux points from a higher-pressure
% donor to a lower-pressure receiver.  Sorting by pressure is therefore
% valid.  In floating-point arithmetic, however, two theoretically equal
% pressures can differ by roundoff, so a fixed tie-breaking order may place
% a receiver before its donor.  The one-dimensional directed flux graph is
% a directed forest and hence acyclic.  We therefore build its topological
% order directly from the signs of the computed fluxes.

N = numel(nnew);
lambda = dt/h;

if size(nhatSpecies,1) ~= N || size(nhatSpecies,2) ~= 2
    error('nhatSpecies must be N-by-2.');
end
if numel(J) ~= N+1
    error('J must have N+1 interface values.');
end

J = J(:);

Jplus = max(J,0);
Jminus = min(J,0);

diagonalMass = nnew ...
    +lambda*(Jplus(2:N+1)-Jminus(1:N));

massScale = max([1;abs(nnew(:));abs(nhatSpecies(:));lambda*abs(J(:))]);

massTol = max( ...
    par.remap_mass_tol*massScale, ...
    100*eps(massScale)*massScale);

negativeTol = max( ...
    par.roundoff_tol*massScale, ...
    100*eps(massScale)*massScale);

if any(diagonalMass < -negativeTol)
    error('Negative remapping diagonal: %.3e.',min(diagonalMass));
end

diagonalMass(diagonalMass < 0) = 0;
active = diagonalMass > massTol;

% Ignore only interface fluxes whose transported mass over one time step
% is below the remapping mass tolerance.  This removes irrelevant
% roundoff-level directions without changing the resolved balance.
Jwork = J;
tinyFlux = lambda*abs(Jwork) <= massTol;
Jwork(tinyFlux) = 0;

%% Build the directed dependency graph
% A receiver depends on every upwind donor carrying mass into it.
indegree = zeros(N,1);

for j = 1:N
    if ~active(j)
        continue;
    end

    % Inflow from the left cell j-1.
    if j > 1 && Jwork(j) > 0
        if ~active(j-1)
            error(['Active cell %d receives a resolved left inflow from ', ...
                'inactive cell %d.'],j,j-1);
        end
        indegree(j) = indegree(j)+1;
    end

    % Inflow from the right cell j+1.
    if j < N && Jwork(j+1) < 0
        if ~active(j+1)
            error(['Active cell %d receives a resolved right inflow from ', ...
                'inactive cell %d.'],j,j+1);
        end
        indegree(j) = indegree(j)+1;
    end
end

queue = zeros(N,1);
seed = find(active & indegree == 0);
tail = numel(seed);
queue(1:tail) = seed;
head = 1;

theta = zeros(N,2);
processed = false(N,1);

maxCompatibilityDefect = 0;
processedCount = 0;

%% Topological sweep
while head <= tail
    j = queue(head);
    head = head+1;

    if processed(j)
        continue;
    end

    rhs = nhatSpecies(j,:);

    if j > 1 && Jwork(j) > 0
        if ~processed(j-1)
            error('Left inflow donor is unavailable at cell %d.',j);
        end
        rhs = rhs+lambda*Jwork(j)*theta(j-1,:);
    end

    if j < N && Jwork(j+1) < 0
        if ~processed(j+1)
            error('Right inflow donor is unavailable at cell %d.',j);
        end
        rhs = rhs-lambda*Jwork(j+1)*theta(j+1,:);
    end

    localScale = max([massTol,abs(rhs),diagonalMass(j)]);
    localNegativeTol = max( ...
        par.roundoff_tol*localScale, ...
        100*eps(localScale)*localScale);

    if any(rhs < -localNegativeTol)
        error('Negative remapping right-hand side at cell %d.',j);
    end

    rhs(rhs < 0) = 0;
    rhsMass = sum(rhs);

    maxCompatibilityDefect = max( ...
        maxCompatibilityDefect, ...
        abs(diagonalMass(j)-rhsMass));

    if rhsMass > massTol
        theta(j,:) = rhs/rhsMass;
    else
        localMass = sum(nhatSpecies(j,:));

        if localMass > massTol
            theta(j,:) = nhatSpecies(j,:)/localMass;
        else
            theta(j,:) = [0.5,0.5];
        end
    end

    theta(j,theta(j,:) < 0) = 0;
    rowSum = sum(theta(j,:));

    if ~(isfinite(rowSum) && rowSum > 0)
        error('Invalid composition at cell %d.',j);
    end

    theta(j,:) = theta(j,:)/rowSum;
    processed(j) = true;
    processedCount = processedCount+1;

    % Outflow from j to j+1.
    if j < N && Jwork(j+1) > 0 && active(j+1)
        indegree(j+1) = indegree(j+1)-1;

        if indegree(j+1) == 0
            tail = tail+1;
            queue(tail) = j+1;
        elseif indegree(j+1) < 0
            error('Negative dependency count at cell %d.',j+1);
        end
    end

    % Outflow from j to j-1.
    if j > 1 && Jwork(j) < 0 && active(j-1)
        indegree(j-1) = indegree(j-1)-1;

        if indegree(j-1) == 0
            tail = tail+1;
            queue(tail) = j-1;
        elseif indegree(j-1) < 0
            error('Negative dependency count at cell %d.',j-1);
        end
    end
end

if processedCount ~= nnz(active)
    unresolved = find(active & ~processed);
    error(['The directed remapping graph was not fully resolved. ', ...
        'First unresolved cell: %d.'],unresolved(1));
end

theta(~active,:) = 0;

activeTheta = theta(active,:);

if isempty(activeTheta)
    simplexError = 0;
    overshoot = 0;
else
    simplexError = max(abs(sum(activeTheta,2)-1));
    overshoot = max([0;-activeTheta(:);activeTheta(:)-1]);
end

if simplexError > par.remap_bound_tol ...
        || overshoot > par.remap_bound_tol
    error('Remapping violates the simplex.');
end

n1new = nnew.*theta(:,1);
n2new = nnew-n1new;
species = [n1new,n2new];

Jspecies = zeros(N+1,2);

for e = 2:N
    if Jwork(e) > 0
        Jspecies(e,:) = Jwork(e)*theta(e-1,:);
    elseif Jwork(e) < 0
        Jspecies(e,:) = Jwork(e)*theta(e,:);
    end
end

residual = species ...
    +lambda*(Jspecies(2:N+1,:)-Jspecies(1:N,:)) ...
    -nhatSpecies;

reconstructionError = max(abs(sum(species,2)-nnew));
speciesResidual = max(abs(residual(:)));

residualScale = 1+max(abs(nhatSpecies(:)));

if reconstructionError > par.reconstruction_tol
    error('Reconstruction error %.3e.',reconstructionError);
end

if speciesResidual > par.remap_residual_tol*residualScale
    message = sprintf( ...
        ['Species remapping residual %.3e; compatibility defect %.3e.'], ...
        speciesResidual,maxCompatibilityDefect);

    if par.abort_on_remap_failure
        error(message);
    else
        warning('%s',message);
    end
end

info.simplex_error = simplexError;
info.reconstruction_error = reconstructionError;
info.species_residual = speciesResidual;
info.compatibility_defect = maxCompatibilityDefect;
info.ignored_flux_mass = max(lambda*abs(J-Jwork));
end

%% ========================================================================
% Diagnostics and exact solution
% ========================================================================

function R = exact_outer_boundary(t,R0)
R = asinh(exp(t).*sinh(R0));
end

function X = exact_internal_interface(R,R0,X0)
argument = tanh(X0/2).*tanh(R/2)./tanh(R0/2);
argument = min(max(argument,-1+10*eps),1-10*eps);
X = 2*atanh(argument);
end

function value = compute_interior_pressure_residual( ...
    p,pInfinity,h,cutoff)

N = numel(p);

d2p = zeros(N,1);
d2p(2:N-1) = ...
    (p(3:N)-2*p(2:N-1)+p(1:N-2))/(h*h);

mask = false(N,1);
mask(2:N-1) = pInfinity(2:N-1) >= cutoff;

value = h*sum(abs(d2p(mask)+1-p(mask)));
end

function value = discrete_l1(v,h)
value = h*sum(abs(v));
end

function y = roundoff_project_interval(y,lo,hi,tol,label)

lowerRoundoff = y < lo & y >= lo-tol;
upperRoundoff = y > hi & y <= hi+tol;

y(lowerRoundoff) = lo;
y(upperRoundoff) = hi;

if any(y < lo | y > hi)
    warning('%s lies outside [%g,%g].',label,lo,hi);
end
end
