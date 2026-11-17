function output_files=reproduce_figure3(data_root,out_dir)
%REPRODUCE_FIGURE3 Redraw final paper Figure 3 without a simulation.
% data_root contains Figure1, Figure2, ...; all outputs go to out_dir.
output_files=utils.render_saved_panels(data_root,out_dir,3);
end
