function [a,p,J,info] = dvd_scalar_corrector(b,sP,dt,par,h)
% Residual r=a-b-dt*L*p+dt*M(p)*(p-pH).
% Jacobian I-dt*L*diag(dp)+dt*diag((M+Mprime*mu)*dp).
N=numel(b); e=.5*(sum(sP(1:end-1,:),2)+sum(sP(2:end,:),2));
L=sparse([1:N-1,2:N,1:N-1,2:N],[2:N,1:N-1,1:N-1,2:N],[e;e;-e;-e]/h^2,N,N);
a=b; info.history=[]; info.converged=false; info.reason='iteration limit';
for it=1:par.newton_maxit
 [r,A,p]=evaluate(a); info.history(end+1)=norm(r,inf);
 if norm(r,inf)<par.newton_tol*(1+norm(b,inf)), info.converged=true; info.reason='root'; break; end
 if strcmp(par.mode,'B')
  grad=A'*r; pg=a-min(max(a-grad,0),par.nH);
  if norm(pg,inf)<par.newton_tol, info.reason='constrained stationary point'; break; end
 end
 d=-(A\r); alpha=1; accepted=false;
 for ls=1:45
  trial=a+alpha*d;
  if strcmp(par.mode,'B'), trial=min(max(trial,0),par.nH); end
  if all(trial>=0) || (strcmp(par.mode,'A') && par.gamma==round(par.gamma))
   rt=evaluate(trial);
   if norm(rt,2)<(1-1e-4*alpha)*norm(r,2), accepted=true; break; end
  end
  alpha=alpha/2;
 end
 if ~accepted && strcmp(par.mode,'B')
  d=-grad/max(1,norm(A,1)*norm(A,inf));alpha=1;
  for ls=1:45
   trial=min(max(a+alpha*d,0),par.nH);rt=evaluate(trial);
   if norm(rt,2)<norm(r,2),accepted=true;break;end
   alpha=alpha/2;
  end
 end
 if ~accepted
  info.reason='domain/line search failure'; break
 end
 a=trial;
end
[r,A,p,M]=evaluate(a); info.residual=norm(r,inf); info.iter=it-1;
grad=A'*r; info.kkt=norm(a-min(max(a-grad,0),par.nH),inf);
info.active=nnz(a<=1e-12|a>=par.nH-1e-12); info.M=M; info.balance=r;
J=[0;-e.*diff(p)/h;0];
 function [r,A,p,M]=evaluate(a)
  [r,A,p,M]=src.dvd_scalar_equation(a,b,sP,dt,par,h,L);
 end
end
