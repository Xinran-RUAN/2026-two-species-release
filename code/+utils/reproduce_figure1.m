function output_files=reproduce_figure1(data_root,out_dir)
%REPRODUCE_FIGURE1 Redraw final paper Figure 1 without a simulation.
% data_root contains Figure1, Figure2, ...; all outputs go to out_dir.
output_files=utils.render_saved_panels(data_root,out_dir,1);
end
