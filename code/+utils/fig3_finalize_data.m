function checks=fig3_finalize_data(out_dir,template_file)
%FIG3_FINALIZE_DATA Build paper plotting arrays and quantify archive differences.
if nargin<2
 paths=utils.paths;template_file=fullfile(paths.data,'Figure3','archived_plot_data.mat');
end
v=load(template_file,'plot_data');plot_data=v.plot_data;archived=plot_data;
a=load(fullfile(out_dir,'production_history_m2.mat'),'result');two=a.result;
a=load(fullfile(out_dir,'production_history_m3.mat'),'result');three=a.result;
ix2=find(abs(two.history.time-1)<1e-12,1);ix3=find(abs(three.history.time-1)<1e-12,1);
U2=two.history.U(:,:,ix2);U3=three.history.U(:,:,ix3);
comparison=struct([]);
for p=1:3
 for j=1:numel(plot_data.panels(p).curves)
  c=plot_data.panels(p).curves(j);if ~strcmp(c.role,'data'),continue;end
  name=c.name;
  if p<=2
   U=U2;if p==2,U=U3;end
   if strcmp(name,'$n$'),y=sum(U,2);id='total_density';else
    si=sscanf(strrep(name,'$',''),'n_%d');y=U(:,si);id=sprintf('species_%d',si);
   end
  else
   colors=[0 .447 .741;.85 .325 .098;0 .55 1;1 .12 .12;.2 .62 .2];
   index=find(vecnorm(colors-c.color,2,2)<1e-6,1);assert(~isempty(index));
   if index<=2,y=two.diagnostics.mass_fraction(:,index);id=sprintf('two_species_fraction_%d',index);
   else,y=three.diagnostics.mass_fraction(:,index-2);id=sprintf('three_species_fraction_%d',index-2);end
  end
  assert(numel(y)==numel(c.y));
  diff=y(:)-c.y(:);
  row=struct('panel',p,'curve_index',j,'curve_id',id,'max_abs_difference',max(abs(diff)), ...
    'mean_abs_difference',mean(abs(diff)),'relative_l2_difference',norm(diff)/max(norm(c.y(:)),realmin));
  if isempty(comparison),comparison=row;else,comparison(end+1)=row;end
  c.y=reshape(y,size(c.y));plot_data.panels(p).curves(j)=c;
  writematrix([c.x(:),c.y(:)],fullfile(out_dir,sprintf('panel%02d_curve%02d.csv',p,j)));
 end
end
plot_data.provenance='Recomputed independently for two and three species using the paper parameters and Newton tolerance 1e-12.';
save(fullfile(out_dir,'plot_data.mat'),'plot_data','-v7');
writetable(struct2table(comparison),fullfile(out_dir,'production_archive_comparison.csv'));
aggregation_density=max(abs(two.history.n(:)-three.history.n(:)));
aggregation_pressure=max(abs(two.history.p(:)-three.history.p(:)));
A=three.history.U(:,1,:)+three.history.U(:,2,:);B=two.history.U(:,1,:);
aggregation_fast=max(abs(A(:)-B(:)));
A=three.history.U(:,3,:);B=two.history.U(:,2,:);aggregation_slow=max(abs(A(:)-B(:)));
checks=struct('two_species',summary(two),'three_species',summary(three), ...
 'max_archived_curve_difference',max([comparison.max_abs_difference]), ...
 'aggregation_density_max_abs_difference',aggregation_density, ...
 'aggregation_pressure_max_abs_difference',aggregation_pressure, ...
 'aggregation_fast_species_max_abs_difference',aggregation_fast, ...
 'aggregation_slow_species_max_abs_difference',aggregation_slow, ...
 'newton_tolerance',two.parameters.newton_tol,'archived_two_species_newton_tolerance',1e-11);
save(fullfile(out_dir,'production_checks.mat'),'checks','-v7');
S=[struct2table(checks.two_species);struct2table(checks.three_species)];S.species_count=[2;3];
writetable(S,fullfile(out_dir,'production_structure_summary.csv'));
fprintf('%s\n',jsonencode(checks,'PrettyPrint',true));
disp(struct2table(comparison));disp(checks);
end
function s=summary(r)
d=r.diagnostics;
s=struct('steps',numel(d.time)-1,'snapshots',numel(r.history.time), ...
 'min_species',min(d.min_species),'max_pressure',max(d.max_pressure), ...
 'max_reconstruction',max(d.reconstruction),'max_simplex',max(d.simplex), ...
 'max_newton_residual',max(d.newton_residual),'max_scalar_balance',max(d.scalar_balance), ...
 'max_species_balance',max(d.species_balance),'max_mass_balance',max(d.mass_balance), ...
 'max_allocation',max(d.allocation),'max_positive_energy_jump',max([0;diff(d.energy)]), ...
 'max_energy_defect',max(d.energy_defect),'max_species_mass_decrease',max([0;reshape(-diff(d.mass_species,1,1),[],1)]), ...
 'max_total_mass_decrease',max([0;-diff(d.mass)]),'wall_seconds',r.wall_seconds);
assert(s.min_species>=0 && s.max_pressure<=1+1e-12 && s.max_positive_energy_jump<=1e-12);
assert(s.max_total_mass_decrease<=1e-12 && s.max_species_mass_decrease<=1e-12);
end
