function finalize_accuracy_case(out_dir,template_file)
%FINALIZE_ACCURACY_CASE Convert recomputed errors to the paper plotting format.
v=load(template_file,'plot_data');plot_data=v.plot_data;
T=readtable(fullfile(out_dir,'temporal_errors.csv'));S=readtable(fullfile(out_dir,'spatial_errors.csv'));
for p=1:numel(plot_data.panels)
 for j=1:numel(plot_data.panels(p).curves)
  c=plot_data.panels(p).curves(j);if ~strcmp(c.role,'data'),continue;end
  if contains(lower(c.name),'near-vacuum'),name='vacuum';else,name='smooth';end
  if p==1
   method='first';if contains(c.name,'PC'),method='dvd';end
   rows=T(strcmp(T.case_name,name)&strcmp(T.method,method),:);x=rows.dt;y=rows.Esp;
  else
   rows=S(strcmp(S.case_name,name),:);x=1./rows.N;y=rows.E_sp_space;
  end
  for k=1:numel(c.x)
   ix=find(abs(x-c.x(k))<1e-14,1);assert(~isempty(ix));c.y(k)=y(ix);
  end
  plot_data.panels(p).curves(j)=c;
 end
end
save(fullfile(out_dir,'plot_data.mat'),'plot_data','-v7');
end
