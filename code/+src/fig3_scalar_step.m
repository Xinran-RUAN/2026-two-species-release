% Original scalar density/Newton equations copied without numerical changes.
function [nnew,pnew,J,info] = fig3_scalar_step(nold,Mmob,par,h)
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
