function [s,info,bal]=fallback_onsager_step(old,dt,par,h)
% Full original Onsager step. Existing local Newton/reaction/remap helpers
% are copied verbatim from simple_predictor_step; only outputs/timing added.
par.dt=dt;par=local_defaults(par);par.abort_on_newton_failure=false;
b=sum(old,2);M=-(old(:,1).*par.sigma1(b.^par.gamma)+old(:,2).*par.sigma2(b.^par.gamma));
scale=max([1;abs(M)]);tol=max(par.roundoff_tol*scale,100*eps(scale)*scale);
if min(M)<-tol,error('Negative Onsager reaction mobility');end
mclip=nnz(M<0);M(M<0)=0;
t=tic;
if isfield(par,'onsager_scalar_solver')
 [a,p,J,ni]=par.onsager_scalar_solver(b,M,dt,par,h);
else
 [a,p,J,ni]=density_step(b,M,par,h);
end
scalar_seconds=toc(t);
info=struct('failed',false,'reason','','newton',ni,'scalar_seconds',scalar_seconds,'clip_count',mclip,'initial_projection_count',nnz(b<0|b>par.nH),'dt',dt);
s=[];bal=struct;
if ~ni.converged,info.failed=true;info.reason='Onsager scalar residual not converged';return;end
hat=b-dt*M.*(p-par.pH);
[z,ri]=local_reaction_allocation(old(:,1),old(:,2),hat,p,par);
if isfield(ni,'pressure_order'),par.pressure_order=ni.pressure_order;end
[s,~,ti]=conservative_species_remap(a,z,J,p,dt,h,par);
info.clip_count=info.clip_count+ri.clip_count+ti.clip_count;
qs=-M.*(p-par.pH);q=(z-old)/dt;Js=ti.Jspecies;
rs=s-old+dt/h*diff(Js)-dt*q;rho=a-b+dt/h*diff(J)-dt*qs;
E=@(n)h*sum(n.^(par.gamma+1)/(par.gamma+1)-par.pH*n);
e=.5*(b(1:end-1)+b(2:end));mu=p-par.pH;
D=h*(sum(e.*(diff(p)/h).^2)+sum(M.*mu.^2));
bal=struct('scalar',a,'pressure',p,'J',J,'Js',Js,'q',q,'qs',qs,'rho',rho,'M',M,'dissipation',D,'energy_change',E(sum(s,2))-E(b),'energy_defect',E(sum(s,2))-E(b)+dt*D,'reconstruction',max(abs(sum(s,2)-a)),'reaction_defect',max(abs(sum(q,2)-qs)),'flux_defect',max(abs(sum(Js,2)-J)),'species_balance',max(abs(rs(:))),'mass_balance',abs(h*sum(sum(s,2)-b)-dt*h*sum(qs)));
if any(~isfinite(s(:))) || any(s(:)<0) || any(sum(s,2)>par.nH)
 info.failed=true;info.reason='Full Onsager step returned invalid raw state';
end
end
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

%% ========================================================================
function [nnew,pnew,J,info] = density_step(nold,Mmob,par,h)
N = length(nold);
dt = par.dt;
n = min(max(nold,0),par.nH);

edge_n = zeros(N+1,1);
edge_n(2:N) = 0.5*(nold(1:N-1)+nold(2:N));

info.projection_trial_count=0;info.history=[];
info.converged = false;
info.iter = 0;
info.residual = Inf;
stop_tol = par.newton_tol*(1+norm(nold,inf));

for it = 1:par.newton_maxit
    [R,Jmat] = density_residual_jac(n,nold,edge_n,Mmob,par,h);
    res_inf = norm(R,inf);info.history(end+1)=res_inf;
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
        rawtrial=n+lambda*dn;
        info.projection_trial_count=info.projection_trial_count+nnz(rawtrial<0|rawtrial>par.nH);
        ntrial = min(max(rawtrial,0),par.nH);
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
clip_count=0;

nhat_scale = max([1;abs(nhat(:))]);
nhat_tol = max(par.roundoff_tol*nhat_scale, ...
    100*eps(nhat_scale)*nhat_scale);
if min(nhat) < -nhat_tol
    error('Negative source-updated total density: min(nhat)=%.3e.', ...
        min(nhat));
end
clip_count=clip_count+nnz(nhat<0);
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
    clip_count=clip_count+nnz(k12<0)+nnz(k21<0);
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
    clip_count=clip_count+nnz(y1<0)+nnz(y2<0);
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
info.clip_count=clip_count;info.normalized_cells=nnz(active);
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
clip_count=nnz(d<0);
d(d<0) = 0;
active = d > mass_tol;

if isfield(par,'pressure_order'),order=par.pressure_order;
else,[~,order]=sort(pnew,'descend');end
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
    clip_count=clip_count+nnz(rhs<0);
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
info.clip_count=clip_count;info.normalized_cells=nnz(active);
info.Jspecies=Jspecies;
end
