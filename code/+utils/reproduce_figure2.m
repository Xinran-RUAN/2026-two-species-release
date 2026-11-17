function output_files=reproduce_figure2(data_root,out_dir)
%REPRODUCE_FIGURE2 Redraw final paper Figure 2 without a simulation.
% data_root contains Figure1, Figure2, ...; all outputs go to out_dir.
output_files=utils.render_saved_panels(data_root,out_dir,2);
end
