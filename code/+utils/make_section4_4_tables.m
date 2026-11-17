function make_section4_4_tables(dataFile,outFolder)
%MAKE_SECTION4_4_TABLES
% Direct-run table generator for Section 4.4.
%
% Usage: utils.make_section4_4_tables(checkpoint_file, out_dir).
% An explicit numeric MAT checkpoint is required; no search is performed.
% Outputs are written to out_dir, or to a temporary directory if omitted.
%
% It recomputes and exports
%   e_n^gamma, e_p^gamma, e_R^gamma, e_X^gamma,
%   their convergence rates,
%   C_gamma and the interior pressure defect I_gamma.
%
% Outputs:
%   tables/table2_errors.tex
%   tables/table3_structure.tex
%   tables/table2_errors.csv
%   tables/table3_structure.csv

clc;

%% ========================================================================
% Automatically locate the newest data file
% ========================================================================

if nargin < 1 || isempty(dataFile), error('An explicit numeric data file is required.'); end

fprintf('Using data file:\n  %s\n\n',dataFile);

S = load(dataFile);

required = { ...
    'x','h','gammas','par','reference', ...
    'final_n','final_p','final_n1','final_n2'};

check_fields(S,required,dataFile);

x = S.x(:);
h = S.h;
gammas = S.gammas(:).';

ng = numel(gammas);

if size(S.final_n,1) ~= numel(x) || size(S.final_n,2) ~= ng
    error('final_n has incompatible dimensions.');
end

if size(S.final_p,1) ~= numel(x) || size(S.final_p,2) ~= ng
    error('final_p has incompatible dimensions.');
end

if size(S.final_n1,1) ~= numel(x) || size(S.final_n1,2) ~= ng
    error('final_n1 has incompatible dimensions.');
end

if size(S.final_n2,1) ~= numel(x) || size(S.final_n2,2) ~= ng
    error('final_n2 has incompatible dimensions.');
end

% Use only completed gamma runs when the checkpoint contains unfinished
% columns. If no completed flag is stored, require finite final fields.
if isfield(S,'completed')
    valid = logical(S.completed(:).');
else
    valid = all(isfinite(S.final_n),1) ...
        & all(isfinite(S.final_p),1) ...
        & all(isfinite(S.final_n1),1) ...
        & all(isfinite(S.final_n2),1);
end

if ~any(valid)
    error('No completed gamma run is available in the data file.');
end

if any(~valid)
    fprintf('Ignoring unfinished gamma values:');
    fprintf(' %g',gammas(~valid));
    fprintf('\n\n');
end

gammas = gammas(valid);
finalN = S.final_n(:,valid);
finalP = S.final_p(:,valid);
finalN1 = S.final_n1(:,valid);
finalN2 = S.final_n2(:,valid);

ng = numel(gammas);

%% ========================================================================
% Recompute all diagnostics
% ========================================================================

errorN = NaN(1,ng);
errorP = NaN(1,ng);
errorR = NaN(1,ng);
errorX = NaN(1,ng);

Rnum = NaN(1,ng);
Xnum = NaN(1,ng);

Cgamma = NaN(1,ng);
Hgamma = NaN(1,ng);

for ig = 1:ng
    n = finalN(:,ig);
    p = finalP(:,ig);
    n1 = finalN1(:,ig);
    n2 = finalN2(:,ig);

    errorN(ig) = h*sum(abs(n-S.reference.n_T));
    errorP(ig) = h*sum(abs(p-S.reference.p_T));

    Rnum(ig) = 0.5*h*sum(n);
    Xnum(ig) = 0.5*h*sum(n1-n2);

    errorR(ig) = abs(Rnum(ig)-S.reference.R_T);
    errorX(ig) = abs(Xnum(ig)-S.reference.X_T);

    Cgamma(ig) = h*sum(abs(p.*(1-n)));

    Hgamma(ig) = interior_pressure_residual( ...
        p,S.reference.p_T,h, ...
        S.par.interior_pressure_cutoff);
end

rateN = convergence_rates(errorN,gammas);
rateP = convergence_rates(errorP,gammas);
rateR = convergence_rates(errorR,gammas);
rateX = convergence_rates(errorX,gammas);

%% ========================================================================
% Output folder
% ========================================================================

if nargin < 2 || isempty(outFolder), outFolder=fullfile(tempdir,'two_species_redraw','stiff_tables'); end

if exist(outFolder,'dir') ~= 7
    mkdir(outFolder);
end

table3Tex = fullfile(outFolder,'table2_errors.tex');
table4Tex = fullfile(outFolder,'table3_structure.tex');

table3Csv = fullfile(outFolder,'table2_errors.csv');
table4Csv = fullfile(outFolder,'table3_structure.csv');

%% ========================================================================
% Write LaTeX Table 3
% ========================================================================

fid = fopen(table3Tex,'w');

if fid < 0
    error('Cannot open %s.',table3Tex);
end

fprintf(fid,'\\begin{table}[htbp]\n');
fprintf(fid,'  \\centering\n');
fprintf(fid,['  \\caption{Field errors, interface errors, and convergence ', ...
    'rates with respect to $1/\\gamma$.}\n']);
fprintf(fid,'  \\label{tab:stiff-pressure-errors}\n');
fprintf(fid,'  \\begin{tabular}{c%s}\n',repmat('c',1,ng));
fprintf(fid,'    \\toprule\n');

fprintf(fid,'    $\\gamma$');
for ig = 1:ng
    fprintf(fid,' & %g',gammas(ig));
end
fprintf(fid,' \\\\\n');

fprintf(fid,'    \\midrule\n');

write_error_row(fid,'$e_n^\\gamma$',errorN);
write_rate_row(fid,'rate',rateN);

write_error_row(fid,'$e_p^\\gamma$',errorP);
write_rate_row(fid,'rate',rateP);

write_error_row(fid,'$e_R^\\gamma$',errorR);
write_rate_row(fid,'rate',rateR);

write_error_row(fid,'$e_X^\\gamma$',errorX);
write_rate_row(fid,'rate',rateX);

fprintf(fid,'    \\bottomrule\n');
fprintf(fid,'  \\end{tabular}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);

%% ========================================================================
% Write LaTeX Table 4
% ========================================================================

fid = fopen(table4Tex,'w');

if fid < 0
    error('Cannot open %s.',table4Tex);
end

fprintf(fid,'\\begin{table}[htbp]\n');
fprintf(fid,'  \\centering\n');
fprintf(fid,'  \\caption{Hele--Shaw structure diagnostics.}\n');
fprintf(fid,'  \\label{tab:stiff-pressure-structure}\n');
fprintf(fid,'  \\begin{tabular}{c%s}\n',repmat('c',1,ng));
fprintf(fid,'    \\toprule\n');

fprintf(fid,'    $\\gamma$');
for ig = 1:ng
    fprintf(fid,' & %g',gammas(ig));
end
fprintf(fid,' \\\\\n');

fprintf(fid,'    \\midrule\n');

write_error_row(fid,'$C_\\gamma$',Cgamma);
write_error_row(fid,'$\\mathcal I_\\gamma$',Hgamma);

fprintf(fid,'    \\bottomrule\n');
fprintf(fid,'  \\end{tabular}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);

%% ========================================================================
% Write CSV files
% ========================================================================

write_csv_table3( ...
    table3Csv,gammas, ...
    errorN,rateN,errorP,rateP,errorR,rateR,errorX,rateX);

write_csv_table4( ...
    table4Csv,gammas,Cgamma,Hgamma);

%% ========================================================================
% Print results directly in the command window
% ========================================================================

fprintf('Table 2: field and interface errors\n');
print_numeric_row('gamma',gammas);
print_numeric_row('e_n',errorN);
print_numeric_row('rate_n',rateN);
print_numeric_row('e_p',errorP);
print_numeric_row('rate_p',rateP);
print_numeric_row('e_R',errorR);
print_numeric_row('rate_R',rateR);
print_numeric_row('e_X',errorX);
print_numeric_row('rate_X',rateX);

fprintf('\nTable 3: Hele--Shaw diagnostics\n');
print_numeric_row('C_gamma',Cgamma);
print_numeric_row('I_gamma',Hgamma);

fprintf('\nGenerated files:\n');
fprintf('  %s\n',table3Tex);
fprintf('  %s\n',table4Tex);
fprintf('  %s\n',table3Csv);
fprintf('  %s\n',table4Csv);
fprintf('\nDone.\n');

end

%% ========================================================================
% Utilities
% ========================================================================

function dataFile=locate_latest_data_file
p=utils.paths;dataFile=fullfile(p.data,'Figure5','stiff_pressure_data.mat');
assert(exist(dataFile,'file')==2,'Saved numerical data were not found.');
end

function check_fields(S,fields,sourceName)

for k = 1:numel(fields)
    if ~isfield(S,fields{k})
        error('Missing field %s in %s.',fields{k},sourceName);
    end
end
end

function rates = convergence_rates(errors,gammas)

rates = NaN(size(errors));

for k = 2:numel(errors)
    rates(k) = log(errors(k-1)/errors(k)) ...
        /log(gammas(k)/gammas(k-1));
end
end

function value = interior_pressure_residual( ...
    p,pInfinity,h,cutoff)

N = numel(p);

d2p = zeros(N,1);
d2p(2:N-1) = ...
    (p(3:N)-2*p(2:N-1)+p(1:N-2))/(h*h);

mask = false(N,1);
mask(2:N-1) = pInfinity(2:N-1) >= cutoff;

value = h*sum(abs(d2p(mask)+1-p(mask)));
end

function write_error_row(fid,label,values)

fprintf(fid,'    %s',label);

for k = 1:numel(values)
    fprintf(fid,' & %s',format_scientific(values(k)));
end

fprintf(fid,' \\\\\n');
end

function write_rate_row(fid,label,values)

fprintf(fid,'    %s',label);

for k = 1:numel(values)
    if isnan(values(k))
        fprintf(fid,' & --');
    else
        fprintf(fid,' & %.3f',values(k));
    end
end

fprintf(fid,' \\\\\n');
end

function text = format_scientific(value)

text = sprintf('%.2E',value);
text = strrep(text,'E-0','E-');
text = strrep(text,'E+0','E+');
end

function write_csv_table3( ...
    fileName,gammas,eN,rN,eP,rP,eR,rR,eX,rX)

fid = fopen(fileName,'w');

if fid < 0
    error('Cannot open %s.',fileName);
end

fprintf(fid,'quantity');

for k = 1:numel(gammas)
    fprintf(fid,',gamma_%g',gammas(k));
end

fprintf(fid,'\n');

write_csv_row(fid,'e_n',eN);
write_csv_row(fid,'rate_n',rN);
write_csv_row(fid,'e_p',eP);
write_csv_row(fid,'rate_p',rP);
write_csv_row(fid,'e_R',eR);
write_csv_row(fid,'rate_R',rR);
write_csv_row(fid,'e_X',eX);
write_csv_row(fid,'rate_X',rX);

fclose(fid);
end

function write_csv_table4(fileName,gammas,Cgamma,Hgamma)

fid = fopen(fileName,'w');

if fid < 0
    error('Cannot open %s.',fileName);
end

fprintf(fid,'quantity');

for k = 1:numel(gammas)
    fprintf(fid,',gamma_%g',gammas(k));
end

fprintf(fid,'\n');

write_csv_row(fid,'C_gamma',Cgamma);
write_csv_row(fid,'I_gamma',Hgamma);

fclose(fid);
end

function write_csv_row(fid,label,values)

fprintf(fid,'%s',label);

for k = 1:numel(values)
    if isnan(values(k))
        fprintf(fid,',');
    else
        fprintf(fid,',%.16e',values(k));
    end
end

fprintf(fid,'\n');
end

function print_numeric_row(label,values)

fprintf('  %-12s',label);

for k = 1:numel(values)
    if isnan(values(k))
        fprintf(' %12s','--');
    else
        fprintf(' %12.4e',values(k));
    end
end

fprintf('\n');
end
