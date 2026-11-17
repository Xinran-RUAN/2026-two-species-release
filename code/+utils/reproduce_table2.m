function reproduce_table2(data_root,out_dir)
%REPRODUCE_TABLE2 Recompute Table 2 from the final numeric fields.
if nargin<2 || isempty(out_dir),out_dir=fullfile(tempdir,'two_species_redraw','Table2');end
if ~exist(out_dir,'dir'),mkdir(out_dir);end
scratch=tempname;mkdir(scratch);cleanup=onCleanup(@()rmdir(scratch,'s'));
utils.make_section4_4_tables(fullfile(data_root,'Table2','stiff_pressure_data.mat'),scratch);
for ext={'csv','tex'}
 copyfile(fullfile(scratch,['table2_errors.' ext{1}]),fullfile(out_dir,['table2_errors.' ext{1}]));
end
end
