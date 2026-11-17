function init=structure_initial
%STRUCTURE_INITIAL Original pointwise initial density and composition.
ell_n=.05;ell_phi=.03;x_interface=.25;
init.n0=@(x).49.*(tanh((x+2)./ell_n)-tanh((x-2)./ell_n));
init.phi1=@(x).5.*(1-tanh((x-x_interface)./ell_phi));
end
