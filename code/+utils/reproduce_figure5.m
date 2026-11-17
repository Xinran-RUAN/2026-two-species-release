function reproduce_figure5(data_root,out_dir)
%REPRODUCE_FIGURE5 Redraw Figure 5 from saved numeric data only.
if nargin<2 || isempty(out_dir),out_dir=fullfile(tempdir,'two_species_redraw','Figure5');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
utils.plot_section4_4_density_pressure(fullfile(data_root,'Figure5','stiff_pressure_data.mat'),out_dir);
end
