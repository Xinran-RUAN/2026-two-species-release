function output_files=reproduce_table1(data_root,out_dir)
%REPRODUCE_TABLE1 Rebuild the efficiency table from recorded timing trials.
% Wall-clock times are historical measurements, not solver-independent
% constants. This function checks their medians against the saved summary.
if ~exist(out_dir,'dir'),mkdir(out_dir);end
D=fullfile(data_root,'Table1');summary=readtable(fullfile(D,'cost.csv'));
rows=cell(4,7);j=0;
for i=1:height(summary)
    g=summary.gamma(i);trials=readtable(fullfile(D,sprintf('cost_g%d.csv',g)));
    for method={'first','dvd'}
        mask=strcmp(trials.method,method{1});r=trials(mask,:);assert(height(r)==3);
        seconds=median(r.seconds);err=r.Esp(1);dt=r.dt(1);
        if strcmp(method{1},'first')
            expected=summary.seconds_first(i);scheme='first-order';
        else
            expected=summary.seconds_dvd(i);scheme='predictor-corrector';
        end
        assert(abs(seconds-expected)<1e-10,'Recorded timing summary differs from trial medians.');
        j=j+1;rows(j,:)={g,scheme,dt,err,seconds,summary.speedup(i),3};
    end
end
T=cell2table(rows,'VariableNames',{'gamma','scheme','time_step','relative_temporal_error','time_seconds','speedup','timing_repetitions'});
csv_path=fullfile(out_dir,'Table1.csv');writetable(T,csv_path);
tex_path=fullfile(out_dir,'Table1.tex');fid=fopen(tex_path,'w');cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%% Regenerated from saved historical timing trials; no simulation was run.\n');
fprintf(fid,'\\begin{tabular}{rlrrrr}\n\\hline\n');
fprintf(fid,'$\\gamma$ & Scheme & $\\tau$ & $\\mathrm{err}_t$ & Time (s) & Speedup \\\\\n\\hline\n');
for i=1:height(T)
    scheme=T.scheme{i};dt=T.time_step(i);e=T.relative_temporal_error(i);s=T.time_seconds(i);
    fprintf(fid,'%d & %s & %.2E & %.2E & %.3g & %.3g \\\\\n',T.gamma(i),scheme,dt,e,s,T.speedup(i));
end
fprintf(fid,'\\hline\n\\end{tabular}\n');
output_files={csv_path,tex_path};
fprintf('Reproduced Table 1 from three recorded timings per scheme and gamma.\n');
end
