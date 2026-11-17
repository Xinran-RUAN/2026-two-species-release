function [n0,phi1_0,phi2_0,n1_0,n2_0]=switching_initial(x,par)
%SWITCHING_INITIAL Original neutral composition and near-vacuum density.
N=par.N;
n0 = 0.49.*( ...
    tanh((x+2.0)./par.ell_n)-tanh((x-2.0)./par.ell_n));
n0 = min(max(n0,0),par.nH);

% Neutral and spatially uniform composition.  Since H_1=H_2, the scalar
% density-pressure evolution is independent of the composition, while the
% pressure-dependent switching terms reshape the hidden species fractions.
phi1_0 = 0.5*ones(N,1);
phi2_0 = 1.0-phi1_0;
n1_0 = phi1_0.*n0;
n2_0 = phi2_0.*n0;

end
