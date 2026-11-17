function [a,w,J,info]=rev_scalar_solve(b,P,dt,p,h,method)
% Increment Newton, with dt-scaled residual and correction criteria.
% p.newton_tol retains its original value 1e-12, now a rate accuracy target.
v=zeros(size(b));[r,A,z]=src.rev_scalar_equation(v,b,P,dt,p,h,method);pre=z.pre;
initial_rho=norm(r,inf);speed=norm(z.rhs,inf);target=dt*p.newton_tol*(1+speed);
expected=dt*speed;state_resolution=eps(max(1,norm(b,inf)));
info=struct('converged',false,'reason','iteration limit','iter',0,'history',[],...
 'trace',struct([]),'projection_trial_count',0,'initial_rho',initial_rho,...
 'expected_increment',expected,'old_speed',speed,'target',target,'roundoff_limited',false);
for it=0:p.newton_maxit
 [r,A,z]=src.rev_scalar_equation(v,b,P,dt,p,h,method,pre);
 ri=norm(r,inf);d=-(A\r);di=norm(d,inf);info.history(end+1)=ri;
 % Residual uncertainty from rounding the increment coordinates themselves.
 % This is the Jacobian applied to their half-ulp grid, not a larger tolerance.
 representation=norm(abs(A)*(0.5*eps(abs(v))),inf);
 residual_budget=target+z.roundoff+representation;
 trace=struct('iteration',it,'rho_inf',ri,'rho_over_dt',ri/dt,'correction',di,...
  'increment',norm(v,inf),'target',target,'roundoff',z.roundoff,'representation',representation,'alpha',0,'line_trials',0);
 first_update_required=it==0 && expected>8*state_resolution;
 correction_target=max(target,(4*p.gamma+16)*eps*norm(v,inf));
 if ~first_update_required && ri<=residual_budget && di<=correction_target
  info.converged=true;info.reason='rate and increment residual';
  info.roundoff_limited=ri>target;info.trace=append(info.trace,trace);break;
 end
 if it==p.newton_maxit || any(~isfinite(d)),info.trace=append(info.trace,trace);break;end
 alpha=1;accepted=false;merit=norm(r,2);
 for ls=1:45
  trial=v+alpha*d;density=b+trial;
  if strcmp(method,'first'),domain=all(density>=0 & density<=p.nH);else,domain=true;end
  if domain
   rt=src.rev_scalar_equation(trial,b,P,dt,p,h,method,pre);
   if all(isfinite(rt)) && norm(rt,2)<(1-1e-4*alpha)*merit,accepted=true;break;end
  end
  alpha=alpha/2;
 end
 trace.alpha=alpha;trace.line_trials=ls;info.trace=append(info.trace,trace);
 if ~accepted,info.reason='increment line search failure';break;end
 v=trial;info.iter=info.iter+1;
end
[r,A,z]=src.rev_scalar_equation(v,b,P,dt,p,h,method,pre);a=b+v;w=z.pressure;J=z.J;
% Distinguish the accurate increment equation and the rounded returned state.
returned=a-b-dt*z.rhs;
info.residual=norm(returned,inf);info.increment_residual=norm(r,inf);
info.rho_over_dt=info.residual/dt;info.increment_rho_over_dt=norm(r,inf)/dt;
info.roundoff=z.roundoff;info.representation=representation;
info.residual_tolerance=target+z.roundoff+representation+state_resolution;
info.increment=norm(v,inf);info.returned_increment=norm(a-b,inf);
info.zero_iteration=info.iter==0;info.below_state_resolution=expected<=8*state_resolution;
info.M=z.M;info.balance=returned;info.gradient=z.gradient;info.mu=z.mu;
info.active=nnz(a<=1e-12|a>=p.nH-1e-12);info.kkt=NaN;
info.final_correction=norm(A\r,inf);
info.increment_vector=v;
% Exact ordering of a=b+v in two floating-point words, for pressure ties.
q=a-b;low=(b-(a-q))+(v-q);
[~,info.pressure_order]=sortrows([-a,-low],[1,2]);
info.state_representation=norm(abs(A)*(0.5*eps(abs(a))),inf);
end
function v=append(v,r)
if isempty(v),v=r;else,v(end+1)=r;end
end
