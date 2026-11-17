function [a,w,J,info]=rev_onsager_scalar(b,M,dt,p,h)
% Reuse the shared increment solve with the original frozen mobility.
q=p;q.sigma1=@(x)-M;q.sigma2=@(x)zeros(size(x));
P=[ones(size(b)),zeros(size(b))];
[a,w,J,info]=src.rev_scalar_solve(b,P,dt,q,h,'first');
end
