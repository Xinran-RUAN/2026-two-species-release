function [Js,valid]=fallback_upwind_flux(U,J)
% First-order only. J is physical flux; J>0 takes the left donor.
N=size(U,1);n=sum(U,2);Js=zeros(N+1,2);valid=true;
if any(~isfinite(J)),valid=false;return;end
for signJ=[1,-1]
 edges=find(signJ*J(2:N)>0)+1;
 donor=edges-1;if signJ<0,donor=edges;end
 if any(n(donor)<=0) || any(~isfinite(n(donor)))
  valid=false;return;
 end
 phi=U(donor,1)./n(donor);
 if any(~isfinite(phi)) || any(phi<0|phi>1),valid=false;return;end
 Js(edges,1)=J(edges).*phi;Js(edges,2)=J(edges)-Js(edges,1);
end
end
