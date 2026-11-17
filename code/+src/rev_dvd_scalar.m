function [a,w,J,info]=rev_dvd_scalar(b,P,dt,p,h)
[a,w,J,info]=src.rev_scalar_solve(b,P,dt,p,h,'dvd');
end
