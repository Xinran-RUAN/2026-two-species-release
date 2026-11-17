function p=accuracy_parameters(kind,g,N,dt,T)
%ACCURACY_PARAMETERS The original experiment parameters with portable solvers.
p=model.fallback_parameters(kind,g,N,dt,T);
p.onsager_scalar_solver=@src.rev_onsager_scalar;
p.dvd_scalar_solver=@src.rev_dvd_scalar;
end
