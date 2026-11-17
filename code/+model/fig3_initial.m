function [U,x,h,rates]=fig3_initial(p,m)
%FIG3_INITIAL Initial fields are sampled at the original cell centers.
h=(p.b-p.a)/p.N;x=p.a+((1:p.N)'-.5)*h;
n=.49*(tanh((x+2)/.05)-tanh((x-2)/.05));
f=.5*(1-tanh((x-.25)/.03));fast=n.*f;
if m==2
 U=[fast,n-fast];rates=[10,2];
elseif m==3
 s=.5*(1+tanh((x+.9)/.05));left=fast.*(1-s);
 % Residual forms are algebraically identical and keep all entries positive.
 U=[left,fast-left,n-fast];rates=[10,10,2];
else,error('This paper experiment uses two or three species.');end
end
