function [p0,n0,n1_0,n2_0]=stiff_initial(x,gamma,par)
%STIFF_INITIAL Original pressure patch and separated populations.
    p0 = max(1-cosh(x)./cosh(par.R0),0);
    n0 = p0.^(1/gamma);

    % Initially separated species.
    n1_0 = n0.*double(x < par.X0);
    n2_0 = n0-n1_0;

end
