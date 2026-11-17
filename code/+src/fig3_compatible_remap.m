% Original implicit pressure-ordered remapping generalized from 2 to m columns.
function [n_species,theta,info] = fig3_compatible_remap( ...
    nnew,nhat_species,J,pnew,dt,h,par)
% Pressure-ordered solution of the implicit conservative upwind remapping.

N = length(nnew);
m = size(nhat_species,2);
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
theta = zeros(N,m);
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
            theta(j,:) = ones(1,m)/m;
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
            theta(j,:) = ones(1,m)/m;
        end
    end
end

% The same reconstruction formula for every species column.
n_species = bsxfun(@times,nnew,theta);

Jspecies = zeros(N+1,m);
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
info.Jspecies = Jspecies;
end

%% ========================================================================
