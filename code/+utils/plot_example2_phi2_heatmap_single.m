function plot_example2_phi2_heatmap_single(dataFile,outFolder)
%PLOT_EXAMPLE2_PHI2_HEATMAP_SINGLE
% Single paper figure for the pressure-induced switching example.
%
% The figure shows the space-time map of the type-2 fraction phi_2.
% Composition is displayed only in the occupied region
%   n(x,t) > vacuumTol.
% The pressure contour p=0.8 identifies the high-pressure core.
%
% Usage: utils.plot_example2_phi2_heatmap_single(checkpoint_file, out_dir)
% An explicit numeric MAT checkpoint is required. No simulation is run.
% Outputs: example2_phi2_heatmap.{png,eps,pdf} under out_dir.

clc; close all;

%% ========================================================================
% User-adjustable plotting parameters
% ========================================================================

% Data and output
opt.outputSubfolder = 'paper_figures';
opt.outputName = 'example2_phi2_heatmap';
opt.outputFormats = {'png','eps','pdf'};
opt.exportResolution = 300;
opt.showFigure = false;
opt.exportFigure = true;

% Figure geometry
opt.figurePosition = [100,100,760,520];
opt.axesPosition = [0.110,0.145,0.700,0.775];
opt.colorbarPosition = [0.850,0.145,0.028,0.775];

% Fonts and axes
opt.fontName = 'Times New Roman';
opt.fontSize = 22;
opt.labelFontSize = 22;
opt.titleFontSize = 20;
opt.legendFontSize = 22;
opt.axesLineWidth = 1.5;
opt.tickInterpreter = 'latex';
opt.textInterpreter = 'latex';
opt.tickDir = 'in';
opt.boxOn = 'on';
opt.gridOn = false;
opt.gridAlpha = 0.12;
opt.legendBox = 'on';

% Domain
opt.xLim = [-8,8];

% Vacuum mask
opt.vacuumTol = [];

% Color scale
% phi_2 is a fraction, so [0,1] is the natural default.
% This can be changed here, but [0,1] is recommended for the paper figure.
opt.phi2ColorLimits = [0,1];
opt.phi2ColorbarTicks = 0:0.2:1;
opt.colormapSize = 256;

% Pressure thresholds
opt.corePressure = [];
opt.frontPressureLow = [];
opt.frontPressureHigh = [];

% Core contour
opt.coreContourColor = [0.00,0.00,0.00];
opt.coreContourLineStyle = '--';
opt.coreContourLineWidth = 2.5;

% Upper boundary of the low-pressure front region
opt.frontHighContourColor = [0.20,0.20,0.20];
opt.frontHighContourLineStyle = '--';
opt.frontHighContourLineWidth = 2.3;

% Lower boundary of the low-pressure front region
opt.frontLowContourColor = [0.20,0.20,0.20];
opt.frontLowContourLineStyle = ':';
opt.frontLowContourLineWidth = 2.3;

% Title and legend
opt.showTitle = true;
opt.titleText = 'space--time map of $\phi_2$';
opt.legendLocation = 'southeast';

%% ========================================================================
% Locate and load data
% ========================================================================

if nargin < 1 || isempty(dataFile)
    dataFile = locate_data_file();
end

if exist(dataFile,'file') ~= 2
    error('Data file does not exist: %s',dataFile);
end

S = load(dataFile);

check_fields(S,{'x'},dataFile);

x = S.x(:);

if isfield(S,'t_hist')
    t = S.t_hist(:);
elseif isfield(S,'snapshot_times')
    t = S.snapshot_times(:);
else
    error('The data file must contain t_hist or snapshot_times.');
end

check_fields(S,{'n_hist','p_hist'},dataFile);

nHist = orient_history(S.n_hist,x,'n_hist');
pHist = orient_history(S.p_hist,x,'p_hist');

if isfield(S,'phi2_hist')
    phi2Hist = orient_history(S.phi2_hist,x,'phi2_hist');
elseif isfield(S,'c_hist')
    phi1Hist = orient_history(S.c_hist,x,'c_hist');
    phi2Hist = 1-phi1Hist;
elseif isfield(S,'phi1_hist')
    phi1Hist = orient_history(S.phi1_hist,x,'phi1_hist');
    phi2Hist = 1-phi1Hist;
else
    error('The data file must contain phi2_hist, c_hist, or phi1_hist.');
end

if size(nHist,2) ~= numel(t) || ...
        size(pHist,2) ~= numel(t) || ...
        size(phi2Hist,2) ~= numel(t)
    error('The saved histories do not match the time vector.');
end

% Read thresholds from the data file when available.
if isempty(opt.vacuumTol)
    opt.vacuumTol = get_parameter(S,'vacuum_tol',1.0e-4);
end
if isempty(opt.corePressure)
    opt.corePressure = get_parameter(S,'core_pressure_threshold',0.8);
end
if isempty(opt.frontPressureLow)
    opt.frontPressureLow = get_parameter(S,'front_pressure_low',0.05);
end
if isempty(opt.frontPressureHigh)
    opt.frontPressureHigh = get_parameter(S,'front_pressure_high',0.30);
end

%% ========================================================================
% Prepare heat map
% ========================================================================

phi2Map = phi2Hist;

% Composition is undefined in vacuum.
phi2Map(nHist <= opt.vacuumTol) = NaN;

%% ========================================================================
% Create figure
% ========================================================================

if opt.showFigure
    visibility = 'on';
else
    visibility = 'off';
end

fig = figure( ...
    'Color','w', ...
    'Visible',visibility, ...
    'Units','pixels', ...
    'Position',opt.figurePosition, ...
    'Renderer','painters');

ax = axes( ...
    'Parent',fig, ...
    'Position',opt.axesPosition, ...
    'Color','w');

disable_axes_toolbar(ax);

img = imagesc(ax,x,t,phi2Map.');
set(ax,'YDir','normal');
set(img,'AlphaData',double(isfinite(phi2Map.')));

hold(ax,'on');

colormap(ax,make_centered_phi2_colormap(opt.colormapSize));
caxis(ax,opt.phi2ColorLimits);

[~,hCore] = contour( ...
    ax,x,t,pHist.', ...
    [opt.corePressure,opt.corePressure], ...
    'Color',opt.coreContourColor, ...
    'LineStyle',opt.coreContourLineStyle, ...
    'LineWidth',opt.coreContourLineWidth);

hold(ax,'off');

xlim(ax,opt.xLim);
ylim(ax,[min(t),max(t)]);

xlabel( ...
    ax,'$x$', ...
    'Interpreter',opt.textInterpreter, ...
    'FontSize',opt.labelFontSize);

ylabel( ...
    ax,'$t$', ...
    'Interpreter',opt.textInterpreter, ...
    'FontSize',opt.labelFontSize);

if opt.showTitle
    title( ...
        ax,opt.titleText, ...
        'Interpreter',opt.textInterpreter, ...
        'FontSize',opt.titleFontSize, ...
        'FontWeight','normal');
end

format_axes(ax,opt);

lgd = legend( ...
    ax,hCore, ...
    {sprintf('$p=%.2g$',opt.corePressure)}, ...
    'Location',opt.legendLocation);

set( ...
    lgd, ...
    'Box',opt.legendBox, ...
    'Interpreter',opt.textInterpreter, ...
    'FontName',opt.fontName, ...
    'FontSize',opt.legendFontSize);

cb = colorbar(ax);

set( ...
    cb, ...
    'Ticks',opt.phi2ColorbarTicks, ...
    'TickLabelInterpreter',opt.tickInterpreter, ...
    'FontSize',opt.fontSize);

ylabel( ...
    cb,'$\phi_2$', ...
    'Interpreter',opt.textInterpreter, ...
    'FontSize',opt.labelFontSize);

% Adding a colorbar changes the axes position.
set(ax,'Position',opt.axesPosition);
set(cb,'Position',opt.colorbarPosition);

%% ========================================================================
% Export
% ========================================================================

if nargin < 2 || isempty(outFolder), outFolder=fullfile(tempdir,'two_species_redraw','Figure4'); end

if exist(outFolder,'dir') ~= 7
    mkdir(outFolder);
end

generatedFiles = {};

if opt.exportFigure
    generatedFiles = save_figure_formats( ...
        fig,outFolder,opt.outputName,opt.outputFormats,opt.exportResolution);
end

%% ========================================================================
% Summary
% ========================================================================

fprintf('\nType-2 space-time heat map\n');
fprintf('  data file: %s\n',dataFile);
fprintf('  vacuum mask: n <= %.4e\n',opt.vacuumTol);
fprintf('  color limits: [%.4g, %.4g]\n', ...
    opt.phi2ColorLimits(1),opt.phi2ColorLimits(2));

if ~isempty(generatedFiles)
    fprintf('  generated files:\n');
    for k = 1:numel(generatedFiles)
        fprintf('    %s\n',generatedFiles{k});
    end
end

fprintf('\n');

if opt.showFigure
    drawnow;
else
    close(fig);
end
end

%% ========================================================================
% Local utilities
% ========================================================================

function dataFile=locate_data_file
p=utils.paths;dataFile=fullfile(p.data,'Figure4','switching_data.mat');
assert(exist(dataFile,'file')==2,'Saved numerical data were not found.');
end

function A = orient_history(A,x,name)

nx = numel(x);

if ~isnumeric(A) || ndims(A) ~= 2
    error('%s must be a numeric 2-D array.',name);
end

if size(A,1) == nx
    return;
elseif size(A,2) == nx
    A = A.';
else
    error('%s has size %d-by-%d, but length(x)=%d.', ...
        name,size(A,1),size(A,2),nx);
end
end

function check_fields(S,fields,sourceName)

missing = {};

for k = 1:numel(fields)
    if ~isfield(S,fields{k})
        missing{end+1} = fields{k}; %#ok<AGROW>
    end
end

if ~isempty(missing)
    error('Missing field(s) in %s: %s', ...
        sourceName,strjoin(missing,', '));
end
end

function value = get_parameter(S,name,defaultValue)

value = defaultValue;

if isfield(S,'par') && isfield(S.par,name)
    value = S.par.(name);
end
end

function disable_axes_toolbar(ax)

try
    axtoolbar(ax,{});
catch
end

try
    disableDefaultInteractivity(ax);
catch
end
end

function format_axes(ax,opt)

set( ...
    ax, ...
    'FontName',opt.fontName, ...
    'FontSize',opt.fontSize, ...
    'LineWidth',opt.axesLineWidth, ...
    'TickLabelInterpreter',opt.tickInterpreter, ...
    'Box',opt.boxOn, ...
    'TickDir',opt.tickDir, ...
    'Color','w', ...
    'Layer','top');

if opt.gridOn
    grid(ax,'on');
    ax.GridAlpha = opt.gridAlpha;
else
    grid(ax,'off');
end
end

function cmap = make_centered_phi2_colormap(ncol)
% Blue below 0.5, nearly white at 0.5, yellow-red above 0.5.

positions = [0.00,0.25,0.50,0.75,1.00];

colors = [ ...
    0.10,0.35,0.75; ...
    0.55,0.80,0.95; ...
    0.96,0.96,0.96; ...
    1.00,0.78,0.25; ...
    0.80,0.05,0.05];

query = linspace(0,1,ncol);

cmap = interp1(positions,colors,query,'linear');
cmap = min(max(cmap,0),1);
end

function files = save_figure_formats( ...
    fig,outFolder,baseName,formats,resolution)

if ischar(formats)
    formats = {formats};
end

set(fig,'Color','w','PaperPositionMode','auto');

files = {};

for k = 1:numel(formats)
    fmt = lower(formats{k});
    fileName = fullfile(outFolder,[baseName,'.',fmt]);

    switch fmt
        case 'fig'
            if exist('savefig','file') == 2
                savefig(fig,fileName);
            else
                hgsave(fig,fileName);
            end

        case 'png'
            print(fig,fileName,'-dpng',sprintf('-r%d',resolution));

        case 'eps'
            set(fig,'Renderer','painters');
            print(fig,fileName,'-depsc','-painters');

        case 'pdf'
            set(fig,'Renderer','painters');
            print(fig,fileName,'-dpdf','-painters','-bestfit');

        otherwise
            error('Unsupported format: %s',formats{k});
    end

    files{end+1,1} = fileName; %#ok<AGROW>
end
end
