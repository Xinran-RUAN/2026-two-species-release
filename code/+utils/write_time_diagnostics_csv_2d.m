function write_time_diagnostics_csv_2d(filename,time,energy,max_p, ...
    reconstruction,simplex,newton_it,newton_res,energy_defect, ...
    scalar_balance,remap_residual)
fid = fopen(filename,'w');
if fid<0
    error('Cannot open diagnostics CSV: %s',filename);
end
fprintf(fid,['time,energy,max_p,reconstruction_error,simplex_error,', ...
    'newton_iterations,newton_residual,energy_inequality_defect,', ...
    'scalar_balance_defect,remap_residual\n']);

nt = numel(time)-1;
for k = 0:nt
    if k==0
        nit = NaN;
        nres = NaN;
        edef = NaN;
        sbal = NaN;
        rres = NaN;
    else
        nit = newton_it(k);
        nres = newton_res(k);
        edef = energy_defect(k);
        sbal = scalar_balance(k);
        rres = remap_residual(k);
    end
    fprintf(fid,['%.16e,%.16e,%.16e,%.16e,%.16e,', ...
        '%.16e,%.16e,%.16e,%.16e,%.16e\n'], ...
        time(k+1),energy(k+1),max_p(k+1), ...
        reconstruction(k+1),simplex(k+1), ...
        nit,nres,edef,sbal,rres);
end
fclose(fid);
end

%% ========================================================================
