function [p,init]=dvd_parameters(kind)
p.a=-1;p.b=1;p.N=200;p.gamma=3;p.pH=1;p.nH=1;p.T=1;p.dt=.001;
p.H1=@(z)2*(1-z).*(1+.2*z);p.H2=@(z)(1-z).*(1+.5*z);
p.sigma1=@(z)-2-.4*z;p.sigma2=@(z)-1-.5*z;
p.dsigma1=@(z)-.4+0*z;p.dsigma2=@(z)-.5+0*z;
p.k12=@(z).2*z;p.k21=@(z).1*(1-z).^2;
p.newton_tol=1e-12;p.newton_maxit=60;p.mode='A';p.progress_every=0;
p.store_solution_history=false;
init.n0=@(x).55+.15*cos(pi*x);init.phi1=@(x).5+.2*cos(2*pi*x);
if strcmp(kind,'structure')
 p.a=-10;p.b=10;p.N=2000;p.gamma=20;p.T=2;p.dt=1e-4;
 p.H1=@(z)10*(1-z);p.H2=@(z)2*(1-z);
 p.sigma1=@(z)-10+0*z;p.sigma2=@(z)-2+0*z;
 p.dsigma1=@(z)0*z;p.dsigma2=@(z)0*z;p.k12=@(z)0*z;p.k21=@(z)0*z;
 init.n0=@(x).49*(tanh((x+2)/.05)-tanh((x-2)/.05));
 init.phi1=@(x).5*(1-tanh((x-.25)/.03));
end
end
