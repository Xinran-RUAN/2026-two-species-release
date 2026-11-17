function [r,A,z]=rev_scalar_equation(v,b,P,dt,p,h,method,pre)
% SAME density equations in increment coordinates; no clipping/projection.
% Precomputation contains only fixed old/predictor values.
g=p.gamma;N=numel(b);a=b+v;
assert(g==round(g) && g>=1,'This verified backend supports positive integer gamma');
if nargin<8
 pre.mu0=(b-p.nH).*src.rev_power_sum(b,p.nH+zeros(size(b)),g);
 pre.grad0=diff(b).*src.rev_power_sum(b(2:end),b(1:end-1),g)/h;
 if strcmp(method,'first')
  pre.e=.5*(b(1:end-1)+b(2:end));w=b.^g;
  pre.M=-(P(:,1).*p.sigma1(w)+P(:,2).*p.sigma2(w));
 else,pre.e=.5*(sum(P(1:end-1,:),2)+sum(P(2:end,:),2));end
 e=pre.e;pre.L=sparse([1:N-1,2:N,1:N-1,2:N],[2:N,1:N-1,1:N-1,2:N],[e;e;-e;-e]/h^2,N,N);
end
if strcmp(method,'first')
 pressure_increment=v.*src.rev_power_sum(a,b,g);dp=g*a.^(g-1);
else
 s=ones(size(a));t=g+zeros(size(a));
 for q=1:g-1,s=s.*a+(q+1)*b.^q;t=t.*a+(g-q)*b.^q;end
 pressure_increment=v.*s/(g+1);dp=t/(g+1);
end
w=b.^g+pressure_increment;mu=pre.mu0+pressure_increment;
% Factor the interface pressure difference itself, so a nearly constant
% pressure is never recovered by subtracting separately evaluated powers.
da=diff(b)+diff(v);db=diff(b);
if strcmp(method,'first')
 c=src.rev_power_sum(a(2:end),a(1:end-1),g);
 gradient=da.*c/h;gradient_scale=abs(da.*c)/h;
else
 ar=a(2:end);al=a(1:end-1);br=b(2:end);bl=b(1:end-1);
 ca=zeros(N-1,1);cb=ca;s=ca;hl=ones(N-1,1);
 for q=1:g
  s=br.*s+bl.^(q-1);
  ca=ar.*ca+hl;cb=ar.*cb+s;
  hl=al.*hl+bl.^q;
 end
 gradient=(da.*ca+db.*cb)/((g+1)*h);
 gradient_scale=(abs(da.*ca)+abs(db.*cb))/((g+1)*h);
end
J=[0;-pre.e.*gradient;0];
if strcmp(method,'first'),M=pre.M;Mp=zeros(size(b));
else
 M=-(P(:,1).*p.sigma1(w)+P(:,2).*p.sigma2(w));
 Mp=-(P(:,1).*p.dsigma1(w)+P(:,2).*p.dsigma2(w));
end
reaction=-M.*mu;rhs=-diff(J)/h+reaction;r=v-dt*rhs;
A=speye(N)-dt*pre.L*spdiags(dp,0,N,N)+dt*spdiags((M+Mp.*mu).*dp,0,N,N);
% Local forward-operation scale: includes flux subtraction cancellation.
flux_scale=[0;abs(pre.e).*gradient_scale;0];
roundoff=(4*g+16)*eps*(abs(v)+dt/h*(flux_scale(1:end-1)+flux_scale(2:end))+dt*abs(reaction));
z=struct('a',a,'pressure',w,'mu',mu,'M',M,'J',J,'rhs',rhs,'gradient',gradient,...
 'roundoff',norm(roundoff,inf),'pre',pre);
end
