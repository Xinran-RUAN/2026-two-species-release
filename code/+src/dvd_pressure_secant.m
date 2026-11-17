function [p,mu,dp] = dvd_pressure_secant(a,b,g,pH)
% Integral secant, stable for equal/near equal and zero endpoints.
if g~=round(g) && (any(a(:)<0)||any(b(:)<0)), error('Negative noninteger DVD domain'); end
p=zeros(size(a)); dp=p;
if g==round(g)
 for k=0:g
  p=p+a.^(g-k).*b.^k/(g+1);
  if g-k>0, dp=dp+(g-k)*a.^(g-k-1).*b.^k/(g+1); end
 end
else
 hi=max(a,b); lo=min(a,b); r=zeros(size(a)); ix=hi>0; r(ix)=lo(ix)./hi(ix);
 z=1-r; near=z<1e-4; far=~near & ix;
 p(far)=hi(far).^g.*(-expm1((g+1)*log(r(far))))./((g+1)*z(far));
 % Around b: integral Taylor through fourth order.
 d=a-b; near=abs(d)<=1e-4*max(a,b) & b>0;
 p(near)=b(near).^g+g/2*b(near).^(g-1).*d(near)+g*(g-1)/6*b(near).^(g-2).*d(near).^2+g*(g-1)*(g-2)/24*b(near).^(g-3).*d(near).^3;
 dp(near)=g/2*b(near).^(g-1)+g*(g-1)/3*b(near).^(g-2).*d(near)+g*(g-1)*(g-2)/8*b(near).^(g-3).*d(near).^2;
 far=ix & ~near;
 dp(far)=(a(far).^g-p(far))./(a(far)-b(far));
end
mu=p-pH;
end
