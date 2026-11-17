function [r,A,p,M]=dvd_scalar_equation(a,b,sP,dt,par,h,L)
N=numel(b);
if nargin<8
 e=.5*(sum(sP(1:end-1,:),2)+sum(sP(2:end,:),2));
 L=sparse([1:N-1,2:N,1:N-1,2:N],[2:N,1:N-1,1:N-1,2:N],[e;e;-e;-e]/h^2,N,N);
end
[p,mu,dp]=src.dvd_pressure_secant(a,b,par.gamma,par.pH);
M=-(sP(:,1).*par.sigma1(p)+sP(:,2).*par.sigma2(p));
Mp=-(sP(:,1).*par.dsigma1(p)+sP(:,2).*par.dsigma2(p));
r=a-b-dt*L*p+dt*M.*mu;
A=speye(N)-dt*L*spdiags(dp,0,N,N)+dt*spdiags((M+Mp.*mu).*dp,0,N,N);
end
