function e=fallback_error(U,ref,h,g)
v=[U,sum(U,2),sum(U,2).^g];r=[ref,sum(ref,2),sum(ref,2).^g];d=v-r;
e.L1=h*sum(abs(d),1);e.L2=sqrt(h*sum(d.^2,1));
e.rel_species=sum(abs(U(:)-ref(:)))/sum(abs(ref(:)));
end
