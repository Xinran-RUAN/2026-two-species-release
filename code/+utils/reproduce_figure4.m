function reproduce_figure4(data_root,out_dir)
%REPRODUCE_FIGURE4 Redraw Figure 4 from saved numeric data only.
if nargin<2 || isempty(out_dir),out_dir=fullfile(tempdir,'two_species_redraw','Figure4');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
utils.plot_example2_phi2_heatmap_single(fullfile(data_root,'Figure4','switching_data.mat'),out_dir);
end
