function [U,d,bal]=fallback_step(old,dt,p,h,method)
% No timestep retries, candidate clipping, theta, projection or MUSCL.
d=struct('failed',false,'reason','','fallback',false,'accepted',false,'main_dt',dt,'predictor_dt',dt/2,'fallback_dt',NaN,'predictor_seconds',0,'corrector_seconds',0,'fallback_seconds',0,'scalar_seconds',0,'newton_calls',0,'newton_iterations',0,'baseline_clips',0,'baseline_projections',0,'candidate',struct);
U=[];bal=struct;d.scalar_diagnostics=struct([]);
if strcmp(method,'first')
 [U,info,bal]=src.fallback_onsager_step(old,dt,p,h);
 d.scalar_diagnostics=compact_newton(d.scalar_diagnostics,info.newton,dt);
 d.failed=info.failed;d.reason=info.reason;d.scalar_seconds=info.scalar_seconds;
 d.newton_calls=1;d.newton_iterations=info.newton.iter;
 d.baseline_clips=info.clip_count;d.baseline_projections=info.newton.projection_trial_count;
 return;
end
t=tic;tc=[];
try
 [P,pi]=src.fallback_onsager_step(old,dt/2,p,h);
 d.scalar_diagnostics=compact_newton(d.scalar_diagnostics,pi.newton,dt/2);
 d.scalar_seconds=pi.scalar_seconds;d.newton_calls=1;d.newton_iterations=pi.newton.iter;
 d.baseline_clips=pi.clip_count;d.baseline_projections=pi.newton.projection_trial_count;
 if pi.failed,d.reason=['predictor: ' pi.reason];else
  d.predictor_seconds=toc(t);tc=tic;b=sum(old,2);
  if isfield(p,'dvd_scalar_solver'),[a,w,J,ni]=p.dvd_scalar_solver(b,P,dt,p,h);
  else,[a,w,J,ni]=src.dvd_scalar_corrector(b,P,dt,p,h);end
  d.scalar_diagnostics=compact_newton(d.scalar_diagnostics,ni,dt);
  d.scalar_seconds=d.scalar_seconds+toc(tc);d.newton_calls=2;d.newton_iterations=d.newton_iterations+ni.iter;
  d.candidate=struct('scalar',a,'scalar_residual',ni.residual,'newton_reason',ni.reason,'min_scalar',min(a),'max_pressure',max(a.^p.gamma));
  if ~ni.converged || ni.residual>p.newton_tol*(1+norm(b,inf))
   d.reason=['scalar solve: ' ni.reason];
  elseif any(~isfinite([a;w;J;ni.M]))
   d.reason='nonfinite scalar candidate';
  elseif any(a<0)
   d.reason='negative scalar root';
  elseif any(a>p.nH)
   d.reason='scalar upper bound';
  else
   [Js,ok]=src.fallback_upwind_flux(P,J);
   if ~ok,d.reason='undefined/nonfinite upwind donor';else
    k12=p.k12(w);k21=p.k21(w);
    q=[(p.H1(w)-k12).*P(:,1)+k21.*P(:,2),k12.*P(:,1)+(p.H2(w)-k21).*P(:,2)];
    qs=-ni.M.*(w-p.pH);H=old+dt*q-dt/h*diff(Js);nH=sum(H,2);
    rho=a-b+dt/h*diff(J)-dt*qs;rs=H-old+dt/h*diff(Js)-dt*q;
    E=@(n)h*sum(n.^(p.gamma+1)/(p.gamma+1)-p.pH*n);
    mu=w-p.pH;e=.5*(sum(P(1:end-1,:),2)+sum(P(2:end,:),2));
    D=h*(sum(e.*(diff(w)/h).^2)+sum(ni.M.*mu.^2));
    bh=struct('scalar',a,'pressure',w,'J',J,'Js',Js,'q',q,'qs',qs,'rho',rho,'M',ni.M,'dissipation',D,'energy_change',E(nH)-E(b),'dvd_identity',E(a)-E(b)+dt*D,'dvd_state_identity',E(nH)-E(b)+dt*D,'dvd_chain_residual',E(a)-E(b)+dt*D-h*sum(mu.*rho),'reconstruction',max(abs(nH-a)),'reaction_defect',max(abs(sum(q,2)-qs)),'flux_defect',max(abs(sum(Js,2)-J)),'species_balance',max(abs(rs(:))),'mass_balance',abs(h*sum(nH-b)-dt*h*sum(qs)));
    d.candidate.min_species=min(H(:));d.candidate.reconstruction=bh.reconstruction;
    d.candidate.raw_species=H;
    if any(~isfinite(H(:))) || any(~isfinite(q(:)))
     d.reason='nonfinite species candidate';
    elseif any(H(:)<0)
     d.reason='negative species candidate';
    elseif any(nH<0) || any(nH>p.nH)
     d.reason='species sum density bound';
    elseif bh.reconstruction>20*p.newton_tol*(1+norm(b,inf)) || bh.species_balance>1e-11 || bh.reaction_defect>1e-11 || bh.flux_defect>1e-11
     d.reason='candidate algebraic balance failure';
    else,U=H;bal=bh;d.accepted=true;end
   end
  end
  d.corrector_seconds=toc(tc);
 end
catch err
 d.reason=['candidate exception: ' err.message];
 if ~isempty(tc),d.corrector_seconds=toc(tc);else,d.predictor_seconds=toc(t);end
end
if d.predictor_seconds==0,d.predictor_seconds=toc(t);end
if ~d.accepted
 d.fallback=true;d.fallback_dt=dt;tf=tic;
 % Deliberately recompute from OLD. No a,w,J,P or H enters this call.
 try,[U,li,bal]=src.fallback_onsager_step(old,dt,p,h);
  d.scalar_diagnostics=compact_newton(d.scalar_diagnostics,li.newton,dt);
  d.scalar_seconds=d.scalar_seconds+li.scalar_seconds;
  d.newton_calls=d.newton_calls+1;d.newton_iterations=d.newton_iterations+li.newton.iter;
  d.baseline_clips=d.baseline_clips+li.clip_count;d.baseline_projections=d.baseline_projections+li.newton.projection_trial_count;
  if li.failed,d.failed=true;d.reason=[d.reason '; fallback failed: ' li.reason];end
 catch err,d.failed=true;d.reason=[d.reason '; fallback exception: ' err.message];end
 d.fallback_seconds=toc(tf);
end
end
function rows=compact_newton(rows,ni,dt)
if ~isfield(ni,'increment_residual'),return;end
r=struct('dt',dt,'initial_rho',ni.initial_rho,'returned_rho',ni.residual,...
 'increment_rho',ni.increment_residual,'rate_residual',ni.rho_over_dt,...
 'increment_rate_residual',ni.increment_rho_over_dt,'expected_increment',ni.expected_increment,...
 'increment',ni.increment,'final_correction',ni.final_correction,...
 'zero_iteration',ni.zero_iteration,'below_state_resolution',ni.below_state_resolution,...
 'roundoff_limited',ni.roundoff_limited);
if isempty(rows),rows=r;else,rows(end+1)=r;end
end
