function reproduce_figure6(data_root,out_dir)
%REPRODUCE_FIGURE6 Draw available Figure 6 rows at stored times t=0,1,2.5.
if nargin<2 || isempty(out_dir),out_dir=fullfile(tempdir,'two_species_redraw','Figure6');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
for gamma=[20,3]
 if ~exist(fullfile(data_root,'Figure6',sprintf('gamma%d_data.mat',gamma)),'file'),continue;end
 S=load(fullfile(data_root,'Figure6',sprintf('gamma%d_data.mat',gamma)));
 target=[0,1,2.5]; ids=zeros(size(target));
 for k=1:numel(target)
  [dist,ids(k)]=min(abs(S.snaps.t-target(k)));
  assert(dist<1e-10,'Required saved snapshot is missing.');
 end
 S.snaps.t=S.snaps.t(ids);S.snaps.n=S.snaps.n(ids);S.snaps.p=S.snaps.p(ids);S.snaps.n_species=S.snaps.n_species(ids);
 if isfield(S.snaps,'phi'),S.snaps.phi=S.snaps.phi(ids);end
 S.n=S.snaps.n{end};S.p=S.snaps.p{end};S.n_species=S.snaps.n_species{end};
 S.phi=zeros(size(S.n_species));active=S.n>0;
 for k=1:size(S.n_species,3)
  a=S.n_species(:,:,k);b=zeros(size(S.n));b(active)=a(active)./S.n(active);S.phi(:,:,k)=b;
 end
 if isfield(S.snaps,'phi'),S.phi=S.snaps.phi{end};end
 S.params.T=target(end); % Displayed final time, selected from saved data.
 rowdir=fullfile(out_dir,sprintf('gamma%d',gamma));
 utils.plot_example4_2d_paper_gamma20(S,rowdir);
 for ext={'pdf','png','eps'}
  src=fullfile(rowdir,['fig_ex4_species_density_evolution.' ext{1}]);
  copyfile(src,fullfile(out_dir,sprintf('figure6_gamma%d.%s',gamma,ext{1})));
 end
end
end
