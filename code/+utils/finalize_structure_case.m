function finalize_structure_case(out_dir,template_file)
%FINALIZE_STRUCTURE_CASE Convert new structure history to paper plotting arrays.
v=load(template_file,'plot_data');plot_data=v.plot_data;
s=load(fullfile(out_dir,'structure_history.mat'));history=s.history;d=s.diagnostics;
for p=1:numel(plot_data.panels)
 for j=1:numel(plot_data.panels(p).curves)
  c=plot_data.panels(p).curves(j);if ~strcmp(c.role,'data'),continue;end
  if contains(plot_data.panels(p).title,'density'),field='n';
  elseif contains(plot_data.panels(p).title,'pressure'),field='p';else,field='';end
  if ~isempty(field)
   t=sscanf(strrep(c.name,'$',''),'t=%f');ix=find(abs(history.time-t)<1e-12,1);
   assert(~isempty(ix));c.y=reshape(history.(field)(:,ix),size(c.y));
  elseif strcmp(c.name,'modified energy'),c.x=reshape(d.time,size(c.x));c.y=reshape(d.energy,size(c.y));
  elseif strcmp(c.name,'total mass'),c.x=reshape(d.time,size(c.x));c.y=reshape(d.mass,size(c.y));end
  plot_data.panels(p).curves(j)=c;
 end
end
save(fullfile(out_dir,'plot_data.mat'),'plot_data','-v7');
end
