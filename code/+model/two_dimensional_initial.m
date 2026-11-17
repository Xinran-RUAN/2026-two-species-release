function [n0,phi0]=two_dimensional_initial(X,Y,par)
%TWO_DIMENSIONAL_INITIAL Original density and spatial composition weights.
[Ny,Nx]=size(X);
r = sqrt(X.^2+Y.^2);
n0 = par.initial_density_height*0.5.* ...
    (1-tanh((r-par.initial_radius)/par.initial_width));
n0 = min(max(n0,0),par.nH);

w1 = par.fast_amplitude*exp(-((X-par.fast_center(1)).^2+ ...
    (Y-par.fast_center(2)).^2)/(2*par.fast_sigma^2));
w3 = par.slow_amplitude*exp(-((X-par.slow_center(1)).^2+ ...
    (Y-par.slow_center(2)).^2)/(2*par.slow_sigma^2));
w2 = ones(size(X));

weights = zeros(Ny,Nx,par.m);
weights(:,:,1) = w1;
weights(:,:,2) = w2;
weights(:,:,3) = w3;
weight_sum = sum(weights,3);

phi0 = zeros(Ny,Nx,par.m);
for s = 1:par.m
    phi0(:,:,s) = weights(:,:,s)./weight_sum;
end
end
