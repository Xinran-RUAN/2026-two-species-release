function plot_section4_4_density_pressure(dataFile,outFolder)
%PLOT_SECTION4_4_DENSITY_PRESSURE
% Two-panel paper figure for Section 4.4.
%
% Left:
%   final total-density profiles n_{gamma,h}(x,T) and the Hele--Shaw limit.
%
% Right:
%   final pressure profiles p_{gamma,h}(x,T) and the Hele--Shaw limit.
%
% Usage:
%   utils.plot_section4_4_density_pressure(checkpoint_file, out_dir)
% An explicit numeric MAT checkpoint is required. No simulation is run.

clc; close all;

%% ========================================================================
% Plot parameters
% ========================================================================

opt.selectedGammas = [10,20,160];

% Figure geometry
opt.figurePosition = [100,100,1180,500];
opt.axesPositionLeft = [0.080,0.165,0.385,0.745];
opt.axesPositionRight = [0.575,0.165,0.385,0.745];

% Fonts and axes
opt.fontName = 'Times New Roman';
opt.fontSize = 24;
opt.labelFontSize = 22;
opt.titleFontSize = 20;
opt.legendFontSize = 22;
opt.axesLineWidth = 1.5;

% Axis limits
opt.xLim = [-2.3,2.3];
opt.nYLim = [0,1.05];
opt.pYLim = [0,0.8];

% Titles and legends
opt.showTitle = true;
opt.densityTitle = 'total density at $T=1$';
opt.pressureTitle = 'pressure profile at $T=1$';
opt.densityLegendLocation = 'south';
opt.pressureLegendLocation = 'south';

% Each gamma curve has independent style parameters.
opt.gamma10Color = [0.0000,0.4470,0.7410];
opt.gamma10LineStyle = '--';
opt.gamma10LineWidth = 3.;

opt.gamma40Color = [0.8500,0.3250,0.0980];
opt.gamma40LineStyle = '-.';
opt.gamma40LineWidth = 3.;

opt.gamma320Color = [0.4660,0.6740,0.1880];
opt.gamma320LineStyle = ':';
opt.gamma320LineWidth = 3.;

% Hele--Shaw limit
opt.exactColor = [0.2,0.2,0.2];
opt.exactLineStyle = '-';
opt.exactLineWidth = 2.8;

% Output
opt.outputSubfolder = 'paper_figures';
opt.outputName = 'section4_4_density_pressure_profiles';
opt.outputFormats = {'png','eps','pdf'};
opt.exportResolution = 300;
opt.showFigure = false;

%% ========================================================================
% Load data
% ========================================================================

if nargin < 1 || isempty(dataFile)
    dataFile = locate_data_file();
end

if exist(dataFile,'file') ~= 2
    error('Data file does not exist: %s',dataFile);
end

S = load(dataFile);

required = {'x','gammas','final_n','final_p','reference','par'};
check_fields(S,required,dataFile);

x = S.x(:);
gammas = S.gammas(:).';

nAll = S.final_n;
pAll = S.final_p;

if size(nAll,1) ~= numel(x) || size(nAll,2) ~= numel(gammas)
    error('final_n has incompatible dimensions.');
end

if size(pAll,1) ~= numel(x) || size(pAll,2) ~= numel(gammas)
    error('final_p has incompatible dimensions.');
end

nExact = S.reference.n_T(:);
pExact = S.reference.p_T(:);

if numel(nExact) ~= numel(x) || numel(pExact) ~= numel(x)
    error('The reference profiles have incompatible dimensions.');
end

ids = zeros(size(opt.selectedGammas));

for q = 1:numel(opt.selectedGammas)
    id = find(gammas == opt.selectedGammas(q),1);

    if isempty(id)
        error('gamma=%g is not available in the data file.', ...
            opt.selectedGammas(q));
    end

    ids(q) = id;
end

styles = { ...
    opt.gamma10Color,opt.gamma10LineStyle,opt.gamma10LineWidth; ...
    opt.gamma40Color,opt.gamma40LineStyle,opt.gamma40LineWidth; ...
    opt.gamma320Color,opt.gamma320LineStyle,opt.gamma320LineWidth};

if size(styles,1) ~= numel(ids)
    error(['The number of curve-style rows must equal the number ', ...
        'of selected gamma values.']);
end

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

%% ------------------------------------------------------------------------
% Left panel: total density
% -------------------------------------------------------------------------

ax1 = axes( ...
    'Parent',fig, ...
    'Position',opt.axesPositionLeft);

hold(ax1,'on');

densityHandles = gobjects(numel(ids)+1,1);
densityLabels = cell(numel(ids)+1,1);

for q = 1:numel(ids)
    densityHandles(q) = plot( ...
        ax1,x,nAll(:,ids(q)), ...
        'Color',styles{q,1}, ...
        'LineStyle',styles{q,2}, ...
        'LineWidth',styles{q,3});

    densityLabels{q} = sprintf('$\\gamma=%g$',gammas(ids(q)));
end

densityHandles(end) = plot( ...
    ax1,x,nExact, ...
    'Color',opt.exactColor, ...
    'LineStyle',opt.exactLineStyle, ...
    'LineWidth',opt.exactLineWidth);

densityLabels{end} = 'Hele--Shaw limit';

% Draw the largest-gamma curve above the Hele--Shaw curve.
% The legend handles and labels are unchanged.
uistack(densityHandles(end-1),'top');

xlim(ax1,opt.xLim);
ylim(ax1,opt.nYLim);

xlabel(ax1,'$x$', ...
    'Interpreter','latex', ...
    'FontSize',opt.labelFontSize);

ylabel(ax1,'$n_\gamma(x,T)$', ...
    'Interpreter','latex', ...
    'FontSize',opt.labelFontSize);

if opt.showTitle
    title(ax1,opt.densityTitle, ...
        'Interpreter','latex', ...
        'FontSize',opt.titleFontSize, ...
        'FontWeight','normal');
end

format_axes(ax1,opt);

lgd1 = legend( ...
    ax1,densityHandles,densityLabels, ...
    'Location',opt.densityLegendLocation);

set(lgd1, ...
    'Interpreter','latex', ...
    'FontName',opt.fontName, ...
    'FontSize',opt.legendFontSize, ...
    'Box','on');

%% ------------------------------------------------------------------------
% Right panel: pressure
% -------------------------------------------------------------------------

ax2 = axes( ...
    'Parent',fig, ...
    'Position',opt.axesPositionRight);

hold(ax2,'on');

pressureHandles = gobjects(numel(ids)+1,1);
pressureLabels = cell(numel(ids)+1,1);

for q = 1:numel(ids)
    pressureHandles(q) = plot( ...
        ax2,x,pAll(:,ids(q)), ...
        'Color',styles{q,1}, ...
        'LineStyle',styles{q,2}, ...
        'LineWidth',styles{q,3});

    pressureLabels{q} = sprintf('$\\gamma=%g$',gammas(ids(q)));
end

pressureHandles(end) = plot( ...
    ax2,x,pExact, ...
    'Color',opt.exactColor, ...
    'LineStyle',opt.exactLineStyle, ...
    'LineWidth',opt.exactLineWidth);

pressureLabels{end} = 'Hele--Shaw limit';

% Draw the largest-gamma curve above the Hele--Shaw curve.
% The legend handles and labels are unchanged.
uistack(pressureHandles(end-1),'top');

xlim(ax2,opt.xLim);
ylim(ax2,opt.pYLim);

xlabel(ax2,'$x$', ...
    'Interpreter','latex', ...
    'FontSize',opt.labelFontSize);

ylabel(ax2,'$p_\gamma(x,T)$', ...
    'Interpreter','latex', ...
    'FontSize',opt.labelFontSize);

if opt.showTitle
    title(ax2,opt.pressureTitle, ...
        'Interpreter','latex', ...
        'FontSize',opt.titleFontSize, ...
        'FontWeight','normal');
end

format_axes(ax2,opt);

lgd2 = legend( ...
    ax2,pressureHandles,pressureLabels, ...
    'Location',opt.pressureLegendLocation);

set(lgd2, ...
    'Interpreter','latex', ...
    'FontName',opt.fontName, ...
    'FontSize',opt.legendFontSize, ...
    'Box','on');

%% ========================================================================
% Export
% ========================================================================

if nargin < 2 || isempty(outFolder), outFolder=fullfile(tempdir,'two_species_redraw','Figure5'); end

if exist(outFolder,'dir') ~= 7
    mkdir(outFolder);
end

save_figure_formats( ...
    fig,outFolder,opt.outputName,opt.outputFormats,opt.exportResolution);

fprintf('\nSection 4.4 density-pressure figure saved in:\n  %s\n', ...
    outFolder);

if ~opt.showFigure
    close(fig);
end
end

%% ========================================================================
% Utilities
% ========================================================================

function dataFile=locate_data_file
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

function format_axes(ax,opt)

set(ax, ...
    'FontName',opt.fontName, ...
    'FontSize',opt.fontSize, ...
    'LineWidth',opt.axesLineWidth, ...
    'TickLabelInterpreter','latex', ...
    'TickDir','in', ...
    'Box','on', ...
    'Layer','top');

grid(ax,'off');
end

function save_figure_formats( ...
    fig,outFolder,baseName,formats,resolution)

set(fig,'Color','w','PaperPositionMode','auto');

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
            error('Unsupported output format: %s',fmt);
    end
end
end
