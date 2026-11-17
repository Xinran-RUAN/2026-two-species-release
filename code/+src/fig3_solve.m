function result=fig3_solve(p,m)
%FIG3_SOLVE First-order scalar Onsager/reaction/implicit remapping scheme.
% No timestep adaptation, new limiter or modification of the model is used.
[U,x,h,rates]=model.fig3_initial(p,m);nt=round(p.T/p.dt);
assert(abs(nt*p.dt-p.T)<1e-13);time=(0:nt)'*p.dt;
steps=unique([0;round(p.snapshot_times(:)/p.dt);nt]);steps=steps(steps<=nt);
ns=numel(steps);history.time=steps*p.dt;history.U=zeros(p.N,m,ns);
history.n=zeros(p.N,ns);history.p=zeros(p.N,ns);
d.time=time;d.mass_species=zeros(nt+1,m);d.mass=zeros(nt+1,1);d.energy=zeros(nt+1,1);
d.min_species=zeros(nt+1,1);d.max_pressure=zeros(nt+1,1);d.reconstruction=zeros(nt+1,1);
d.newton_residual=zeros(nt,1);d.newton_iterations=zeros(nt,1);d.energy_defect=zeros(nt,1);
d.scalar_balance=zeros(nt,1);d.species_balance=zeros(nt,1);d.allocation=zeros(nt,1);
d.mass_balance=zeros(nt,1);d.simplex=zeros(nt,1);
[d,history]=store(d,history,1,1,U,p,h);sid=2;wall=tic;
for k=1:nt
 old=U;nold=sum(old,2);M=old*rates(:);Eold=d.energy(k);
 [nnew,pnew,J,info]=src.fig3_scalar_step(nold,M,p,h);
 assert(info.converged,'Newton did not satisfy the requested scaled tolerance.');
 reaction=M.*(p.pH-pnew);nhat=nold+p.dt*reaction;
 % Normalized exp(dt*diag(H_i(pnew))) applied at frozen new pressure.
 active=nold>p.occupied_tol;phi=zeros(p.N,m);phi(:,end)=1;
 z=bsxfun(@times,old(active,:),exp(p.dt*bsxfun(@times,p.pH-pnew(active),rates)));
 den=sum(z,2);assert(all(den>0 & isfinite(den)));
 phi(active,:)=bsxfun(@rdivide,z,den);
 allocated=bsxfun(@times,nhat,phi);
 [U,theta,remap]=src.fig3_compatible_remap(nnew,allocated,J,pnew,p.dt,h,p);
 assert(all(U(:)>=0) && all(isfinite(U(:))));
 [d,history]=store(d,history,k+1,0,U,p,h);
 d.reconstruction(k+1)=max(abs(sum(U,2)-nnew));
 d.newton_residual(k)=info.residual;d.newton_iterations(k)=info.iter;
 d.scalar_balance(k)=max(abs(nnew+p.dt/h*diff(J)-nhat));
 d.species_balance(k)=remap.species_residual;d.simplex(k)=remap.simplex_error;
 d.allocation(k)=max(abs(sum(allocated,2)-nhat));
 d.mass_balance(k)=abs(h*sum(sum(U,2)-nold)-p.dt*h*sum(reaction));
 edge=.5*(nold(1:end-1)+nold(2:end));
 diss=p.dt*h*(sum(edge.*(diff(pnew)/h).^2)+sum(M.*(pnew-p.pH).^2));
 d.energy_defect(k)=d.energy(k+1)+diss-Eold;
 if sid<=ns && k==steps(sid)
  history.U(:,:,sid)=U;history.n(:,sid)=sum(U,2);history.p(:,sid)=sum(U,2).^p.gamma;sid=sid+1;
 end
 if p.progress_every>0 && (mod(k,p.progress_every)==0 || k==nt)
  fprintf('Figure3 m=%d step=%d/%d t=%.4f Newton=%d residual=%.3e Edef=%.3e elapsed=%.1fs\n', ...
   m,k,nt,k*p.dt,info.iter,info.residual,d.energy_defect(k),toc(wall));
 end
end
d.mass_fraction=bsxfun(@rdivide,d.mass_species,d.mass);
result=struct('parameters',p,'species_count',m,'growth_coefficients',rates,'x',x,'h',h, ...
 'history',history,'diagnostics',d,'final_U',U,'wall_seconds',toc(wall));
end
function [d,history]=store(d,history,k,sid,U,p,h)
n=sum(U,2);d.mass_species(k,:)=h*sum(U,1);d.mass(k)=h*sum(n);
d.energy(k)=h*sum(n.^(p.gamma+1)/(p.gamma+1)-p.pH*n);
d.min_species(k)=min(U(:));d.max_pressure(k)=max(n.^p.gamma);
if sid>0,history.U(:,:,sid)=U;history.n(:,sid)=n;history.p(:,sid)=n.^p.gamma;end
end
