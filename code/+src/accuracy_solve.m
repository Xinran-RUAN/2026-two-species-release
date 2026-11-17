function out=accuracy_solve(p,method,snapshot_times)
%ACCURACY_SOLVE Portable fixed-step recorder using the original update loop.
% Cache lookup, source-code snapshots and file writes were removed. Numerical
% updates, diagnostic arithmetic and timed loop remain as in rev_solve.m.
if nargin<3,snapshot_times=[];end
key=sprintf('%s_g%d_N%d_dt%.3g',method,p.gamma,p.N,p.dt);
[U,x,h]=model.fallback_initial(p);nt=round(p.T/p.dt);assert(abs(nt*p.dt-p.T)<32*eps(max(1,p.T)));
E=@(n)h*sum(n.^(p.gamma+1)/(p.gamma+1)-p.pH*n);b=sum(U,2);E0=E(b);mass0=h*sum(b);
z=struct('steps',0,'accepted_steps',0,'fallback_steps',0,'first_steps',0,'min_species',min(U(:)),...
 'max_pressure',max(b.^p.gamma),'pressure_excess',0,'energy_growth',0,'mass_decrease',0,...
 'reconstruction',0,'reaction_defect',0,'flux_defect',0,'species_balance',0,'mass_balance',0,...
 'scalar_residual',0,'dvd_identity',0,'dvd_state_identity',0,'dvd_chain_residual',0,'low_energy_violation',0,...
 'seconds',0,'predictor_seconds',0,'corrector_seconds',0,'fallback_seconds',0,'scalar_seconds',0,...
 'newton_calls',0,'newton_iterations',0,'baseline_clips',0,'baseline_projections',0,...
 'dt_min',Inf,'dt_max',0,'max_rate_residual',0,'max_increment_rate_residual',0,...
 'sum_increment_residual',0,'sum_returned_residual',0,'sum_final_correction',0,...
 'zero_newton',0,'false_zero_newton',0,'roundoff_stops',0);
success=true;reason='';reasons=struct('reason',{},'count',{});rejections=struct([]);
snaps=unique(round(snapshot_times/p.dt));snapshot_U={U};snapshot_t=0;trace=[0,E0,mass0];stride=max(1,ceil(nt/1000));
wall=tic;
for k=1:nt
 old=U;[next,d,bal]=src.fallback_step(old,p.dt,p,h,method);
 for name={'predictor_seconds','corrector_seconds','fallback_seconds','scalar_seconds','newton_calls','newton_iterations','baseline_clips','baseline_projections'}
  f=name{1};z.(f)=z.(f)+d.(f);
 end
 for j=1:numel(d.scalar_diagnostics)
  nd=d.scalar_diagnostics(j);z.max_rate_residual=max(z.max_rate_residual,nd.rate_residual);
  z.max_increment_rate_residual=max(z.max_increment_rate_residual,nd.increment_rate_residual);
  z.sum_increment_residual=z.sum_increment_residual+nd.increment_rho;
  z.sum_returned_residual=z.sum_returned_residual+nd.returned_rho;
  z.sum_final_correction=z.sum_final_correction+nd.final_correction;
  z.zero_newton=z.zero_newton+nd.zero_iteration;
  z.false_zero_newton=z.false_zero_newton+(nd.zero_iteration && ~nd.below_state_resolution && nd.expected_increment>0);
  z.roundoff_stops=z.roundoff_stops+nd.roundoff_limited;
 end
 if d.fallback
  ix=find(strcmp({reasons.reason},d.reason),1);
  if isempty(ix),reasons(end+1)=struct('reason',d.reason,'count',1);else,reasons(ix).count=reasons(ix).count+1;end
  r=struct('time',(k-1)*p.dt,'reason',d.reason,'failed',d.failed,'min_scalar',NaN,'max_pressure',NaN,'pressure_excess',NaN,'min_species',NaN,'scalar_residual',NaN,'reconstruction',NaN);
  for f={'min_scalar','max_pressure','min_species','scalar_residual','reconstruction'}
   n=f{1};if isfield(d.candidate,n),r.(n)=d.candidate.(n);end
  end
  r.pressure_excess=max(0,r.max_pressure-p.pH);rejections=append(rejections,r);
 end
 if d.failed,success=false;reason=d.reason;break;end
 z.steps=k;z.accepted_steps=z.accepted_steps+d.accepted;z.fallback_steps=z.fallback_steps+d.fallback;z.first_steps=z.first_steps+strcmp(method,'first');
 z.dt_min=min(z.dt_min,d.main_dt);z.dt_max=max(z.dt_max,d.main_dt);
 for name={'reconstruction','reaction_defect','flux_defect','species_balance','mass_balance'}
  f=name{1};z.(f)=max(z.(f),bal.(f));
 end
 z.scalar_residual=max(z.scalar_residual,norm(bal.rho,inf));
 if d.accepted
  for name={'dvd_identity','dvd_state_identity','dvd_chain_residual'},f=name{1};z.(f)=max(z.(f),abs(bal.(f)));end
 else,z.low_energy_violation=max(z.low_energy_violation,bal.energy_defect);end
 n=sum(next,2);z.min_species=min(z.min_species,min(next(:)));z.max_pressure=max(z.max_pressure,max(n.^p.gamma));
 z.pressure_excess=max(0,z.max_pressure-p.pH);z.energy_growth=max(z.energy_growth,bal.energy_change);
 z.mass_decrease=max(z.mass_decrease,h*sum(sum(old,2)-n));U=next;
 if any(snaps==k),snapshot_U{end+1}=U;snapshot_t(end+1)=k*p.dt;end
 if mod(k,stride)==0 || k==nt,trace(end+1,:)=[k*p.dt,E(n),h*sum(n)];end
 if mod(k,1000)==0,fprintf('%s %d/%d FB%d elapsed %.1fs\n',key,k,nt,z.fallback_steps,toc(wall));end
end
z.seconds=toc(wall);z.seconds_per_step=z.seconds/max(1,z.steps);z.mean_newton=z.newton_iterations/max(1,z.newton_calls);
z.energy_change=E(sum(U,2))-E0;z.mass_change=h*sum(U(:))-mass0;
out=struct('U',U,'n',sum(U,2),'pressure',sum(U,2).^p.gamma,'h',h,'x',x,'stats',z,...
 'success',success,'reason',reason,'time',z.steps*p.dt,'reasons',reasons,'rejections',rejections,...
 'trace',trace,'snapshot_U',{snapshot_U},'snapshot_times',snapshot_t);
end
function v=append(v,r)
if isempty(v),v=r;else,v(end+1)=r;end
end
