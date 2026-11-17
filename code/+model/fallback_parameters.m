function p=fallback_parameters(kind,g,N,dt,T)
% All required functions are in this portable solver directory.
[p,~]=model.dvd_parameters('smooth');
p.kind=kind;p.gamma=g;p.N=N;p.dt=dt;p.T=T;p.nH=p.pH^(1/g);
p.newton_tol=1e-12;p.newton_maxit=60;p.mode='A';
p.occupied_tol=1e-13;p.remap_mass_tol=1e-14;
p.roundoff_tol=1e-12;p.remap_negative_tol=1e-11;
if ~strcmp(kind,'paper_smooth') && ~strcmp(kind,'paper_structure')
 p.H1=@(z)1-z;p.H2=@(z).6*(1-z);
 p.sigma1=@(z)-1+0*z;p.sigma2=@(z)-.6+0*z;
 p.dsigma1=@(z)0*z;p.dsigma2=@(z)0*z;
 p.k12=@(z).2+0*z;p.k21=@(z).1+0*z;
end
if strcmp(kind,'paper_structure')
 [q,~]=model.dvd_parameters('structure');
 names={'a','b','H1','H2','sigma1','sigma2','dsigma1','dsigma2','k12','k21'};
 for j=1:numel(names),p.(names{j})=q.(names{j});end
end
end
