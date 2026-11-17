function r=s41_error(U,V,h,g)
% Density references are cell averages. Pressure is recomputed AFTER restriction.
e=utils.fallback_error(U,V,h,g);ref=[V,sum(V,2),sum(V,2).^g];
r.En=e.L1(3)/(h*sum(abs(ref(:,3))));r.Esp=e.rel_species;
names={'n1','n2','n','p'};
for j=1:4
 r.(['L1_' names{j}])=e.L1(j);r.(['L2_' names{j}])=e.L2(j);
 r.(['relL1_' names{j}])=e.L1(j)/(h*sum(abs(ref(:,j))));
 r.(['relL2_' names{j}])=e.L2(j)/sqrt(h*sum(ref(:,j).^2));
end
end
