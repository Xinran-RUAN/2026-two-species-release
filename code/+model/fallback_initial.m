function [U,x,h]=fallback_initial(p)
h=(p.b-p.a)/p.N;x=p.a+((1:p.N)'-.5)*h;
if isfield(p,'initial_species'),U=p.initial_species;return;end
if strcmp(p.kind,'smooth') || strcmp(p.kind,'paper_smooth')
 if strcmp(p.kind,'smooth'),c=.55;a=.1;r=.1;else,c=.55;a=.15;r=.2;end
 n=c+a*cs(pi*h/2)*cos(pi*x);
 % Exact cell averages of n(x)*(1/2+r*cos(2*pi*x)).
 n1=c/2+a*(.5+r/2)*cs(pi*h/2)*cos(pi*x)+c*r*cs(pi*h)*cos(2*pi*x)+a*r/2*cs(3*pi*h/2)*cos(3*pi*x);
 U=[n1,n-n1];
else
 z=[-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363];
 w=[.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763]/2;
 q=x+h/2*z;
 if strcmp(p.kind,'paper_structure')
  n=.49*(tanh((q+2)/.05)-tanh((q-2)/.05));phi=.5*(1-tanh((q-.25)/.03));
 else
  n=.98*(abs(q)<.65);phi=.5+.15*cos(2*pi*q);
 end
 U=[(n.*phi)*w',(n.*(1-phi))*w'];
end
end
function v=cs(z)
if z==0,v=1;else,v=sin(z)/z;end
end
