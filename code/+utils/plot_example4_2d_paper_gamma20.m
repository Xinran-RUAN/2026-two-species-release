function plot_example4_2d_paper_gamma20(dataFile, figDir)
%PLOT_EXAMPLE4_2D_PAPER  Publication-quality plot for Example 4.
%
% This script only reads existing data. It does not rerun the simulation.
%
% Main paper outputs:
%   fig_ex4_species_density_evolution.pdf/png/eps
%
% The figure shows RGB-composite evolution of reconstructed species densities
% n_i=phi_i*n at three selected times.
%
% Usage:
%   plot_example4_2d_paper
%   plot_example4_2d_paper('example4_2d_multispecies_output/example4_2d_multispecies_data.mat')
%   plot_example4_2d_paper(dataFile, 'paper_figures')

if nargin < 1 || isempty(dataFile)
    paths=utils.paths;dataFile=fullfile(paths.data,'Figure6','gamma20_data.mat');
end
if nargin < 2 || isempty(figDir)
    figDir = fullfile(tempdir,'two_species_plots','Figure6');
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end
if ~isstruct(dataFile) && ~exist(dataFile,'file'), error('Data file not found.'); end

if isstruct(dataFile), S=dataFile; else, S=load(dataFile); end

x = S.x(:)';
y = S.y(:)';
n = S.n;
p = S.p;
phi = S.phi;

if isfield(S, 'ni')
    ni = S.ni;
else
    [Ny, Nx, m] = size(phi); %#ok<ASGLU>
    ni = zeros(size(phi));
    for s = 1:m
        ni(:,:,s) = phi(:,:,s) .* n;
    end
end

if isfield(S, 'params')
    params = S.params;
else
    params = struct();
end
if isfield(S, 'snaps')
    snaps = S.snaps;
else
    snaps = [];
end
[~, ~, m] = size(phi);

% Growth coefficients for Example 4:
% H_i(p)=a_i(1-p), i=1,2,3.
growthCoeff = get_growth_coeff(params, m);

%% ===================== User-adjustable plot settings =====================
% --- time points shown in the paper figure
plotOpt.targetTimes = [0, 1, 2.5];

% --- level set used as the front
plotOpt.frontLevel = 0.5;

% --- near-vacuum threshold for whitening
plotOpt.vacuumLevel = 0.02;

% --- RGB rendering
plotOpt.densityGamma = 0.75;

% --- displayed window
plotOpt.xlimPlot = [-8, 8];
plotOpt.ylimPlot = [-8, 8];

% --- export options
plotOpt.savePdf = true;
plotOpt.savePng = true;
plotOpt.saveEps = true;
plotOpt.pngResolution = 300;

% --- typography
plotOpt.fontName = 'Times New Roman';
plotOpt.fontSizeAxis = 22;
plotOpt.fontSizeLabel = 22;
plotOpt.fontSizeTitle = 20;
plotOpt.fontSizeLegend = 20;

% --- line styles
plotOpt.axisLineWidth = 1.5;
plotOpt.boxOn = 'on';
plotOpt.tickDir = 'in';
plotOpt.gridOn = false;
plotOpt.lineWidthContour = 1.5;

% --- figure size
plotOpt.figWidthCm = 42.0;
plotOpt.figHeightCm = 17.0;

% --- manual axes layout (normalized units)
plotOpt.leftMargin = 0.04;
plotOpt.rightMargin = 0.12;
plotOpt.bottomMargin = 0.16;
plotOpt.topMargin = 0.06;
plotOpt.hGap = 0.055;

% --- species color key
plotOpt.speciesLegendPosition = [0.92, 0.31, 0.015, 0.44];
plotOpt.speciesLegendFontSize = 22;

plotOpt.useLatex = true;
plotOpt.interpreter = ternary(plotOpt.useLatex, 'latex', 'tex');

if m == 3
    plotOpt.speciesColors = [


    0.80, 0.10, 0.30;   % n1: fast clone, reddish magenta
    0.90, 0.55, 0.15;   % n2: background, orange
    0.00, 0.55, 0.65    % n3: slow clone, teal

    ]; 

    plotOpt.speciesLegendLabels = {'$\phi_1$', '$\phi_2$', '$\phi_3$'};
else
    plotOpt.speciesColors = lines(m);
    plotOpt.speciesLegendLabels = make_species_index_labels(m);
end

snapData = collect_requested_snapshots(snaps, n, phi, params, plotOpt.targetTimes);

%% Paper figure: species density evolution
plot_species_density_evolution(snapData, x, y, plotOpt, figDir);

% Print useful numerical diagnostics for the text.
% Fractions are defined only on occupied cells; vacuum is not a simplex defect.
occupied = n > 0;
simplexErr = abs(sum(phi,3) - 1);
simplexErr = simplexErr(occupied);
if isempty(simplexErr), simplexErr = 0; end
reconErr = abs(sum(ni,3) - n);

fprintf('Paper figure written to: %s\n', figDir);
fprintf('Final diagnostics:\n');
fprintf('  max p                 = %.6g\n', max(p(:)));
fprintf('  min phi               = %.3e\n', min(phi(:)));
fprintf('  active-cell simplex error = %.3e\n', max(simplexErr(:)));
fprintf('  reconstruction error  = %.3e\n', max(reconErr(:)));
fprintf('  growth coefficients   = ');
fprintf('%.4g ', growthCoeff);
fprintf('\n');

end

%% ========================================================================
function plot_species_density_evolution(snapData, x, y, plotOpt, figDir)

nt = numel(snapData.t);

fig = figure('Visible','off','Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, plotOpt.figWidthCm, plotOpt.figHeightCm]);

panelWidth = (1 - plotOpt.leftMargin - plotOpt.rightMargin - 2*plotOpt.hGap) / 3;
panelHeight = 1 - plotOpt.bottomMargin - plotOpt.topMargin;
for q = 1:nt
    offset=1; if isfield(plotOpt,'panelOffset'),offset=plotOpt.panelOffset;end
    pos = [plotOpt.leftMargin + (q-offset)*(panelWidth + plotOpt.hGap), ...
        plotOpt.bottomMargin, panelWidth, panelHeight];
    ax = axes('Parent', fig, 'Position', pos); %#ok<LAXES>
    nSnapshot = snapData.n{q};
    phiSnapshot = snapData.phi{q};
    RGB = make_species_rgb_image(nSnapshot, phiSnapshot, ...
        plotOpt.speciesColors, plotOpt.vacuumLevel, plotOpt.densityGamma);

    image(x, y, RGB);
    set(ax, 'YDir', 'normal');
    axis(ax, 'equal');
    xlim(ax, plotOpt.xlimPlot);
    ylim(ax, plotOpt.ylimPlot);
    axis(ax, 'manual');

    hold(ax, 'on');
    contour(ax, x, y, nSnapshot, [plotOpt.frontLevel plotOpt.frontLevel], ...
        'k-', 'LineWidth', plotOpt.lineWidthContour);
    hold(ax, 'off');

    title(ax, make_time_title(snapData.tTarget(q)), ...
        'Interpreter', plotOpt.interpreter, 'FontSize', plotOpt.fontSizeTitle);
    xlabel(ax, '$x$', 'Interpreter', plotOpt.interpreter, ...
        'FontSize', plotOpt.fontSizeLabel);
    ylabel(ax, '$y$', 'Interpreter', plotOpt.interpreter, ...
        'FontSize', plotOpt.fontSizeLabel);
    format_axes(ax, plotOpt);
    apply_grid(ax, plotOpt);
end

draw_species_color_key(fig, plotOpt);

save_all_formats(fig, fullfile(figDir, 'fig_ex4_species_density_evolution'), ...
    plotOpt.savePdf, plotOpt.savePng, plotOpt.saveEps, plotOpt.pngResolution);
end

function RGB = make_species_rgb_image(nSnap, phiSnap, speciesColors, vacuumLevel, densityGamma)
% Construct a white-background RGB image weighted by species densities.

[Ny, Nx, m] = size(phiSnap);
nShow = min(max(nSnap, 0), 1);
alpha = nShow .^ densityGamma;
RGBcolor = zeros(Ny, Nx, 3);

for s = 1:m
    for ch = 1:3
        RGBcolor(:,:,ch) = RGBcolor(:,:,ch) + phiSnap(:,:,s) .* speciesColors(s,ch);
    end
end

RGB = ones(Ny, Nx, 3);
for ch = 1:3
    RGB(:,:,ch) = (1 - alpha) + alpha .* RGBcolor(:,:,ch);
end

vacuumMask = nSnap < vacuumLevel;
for ch = 1:3
    tmp = RGB(:,:,ch);
    tmp(vacuumMask) = 1;
    RGB(:,:,ch) = tmp;
end

RGB = min(max(RGB, 0), 1);
end

function titleStr = make_time_title(tValue)
if abs(tValue - round(tValue)) < 1e-12
    titleStr = sprintf('$t=%d$', round(tValue));
else
    titleStr = sprintf('$t=%.3g$', tValue);
end
end

function draw_species_color_key(fig, plotOpt)
m = size(plotOpt.speciesColors, 1);
axKey = axes('Parent', fig, 'Position', plotOpt.speciesLegendPosition); %#ok<LAXES>
imagesc(axKey, 1, 1:m, (1:m)');
set(axKey, 'YDir', 'normal', ...
    'XTick', [], ...
    'YTick', 1:m, ...
    'YTickLabel', plotOpt.speciesLegendLabels, ...
    'YAxisLocation', 'right', ...
    'FontName', plotOpt.fontName, ...
    'FontSize', plotOpt.speciesLegendFontSize, ...
    'LineWidth', plotOpt.axisLineWidth, ...
    'Box', 'on', ...
    'Layer', 'top');
try
    set(axKey, 'TickLabelInterpreter', plotOpt.interpreter);
catch
end
xlim(axKey, [0.5, 1.5]);
ylim(axKey, [0.5, m+0.5]);
caxis(axKey, [0.5, m+0.5]);
colormap(axKey, plotOpt.speciesColors);
end

function legendLabels = make_species_index_labels(m)
legendLabels = cell(1,m);
for s = 1:m
    legendLabels{s} = sprintf('$\\phi_%d$', s);
end
end

function snapData = collect_requested_snapshots(snaps, nFinal, phiFinal, params, targetTimes)

if isempty(snaps) || ~isstruct(snaps) || ~isfield(snaps, 'n') || isempty(snaps.n)
    error('snaps.n is required to plot species density evolution.');
end
if ~isfield(snaps, 't') || isempty(snaps.t)
    error('snaps.t is required to select species density snapshots.');
end

tSnap = snaps.t(:);
nCell = snapshot_field_to_cell(snaps.n, size(nFinal), 'snaps.n');
if isfield(snaps, 'phi') && ~isempty(snaps.phi)
    phiCell = snapshot_field_to_cell( ...
        snaps.phi, size(phiFinal(:,:,1)), 'snaps.phi');
elseif isfield(snaps, 'n_species') && ~isempty(snaps.n_species)
    speciesCell = snapshot_field_to_cell( ...
        snaps.n_species, size(phiFinal(:,:,1)), 'snaps.n_species');
    phiCell = cell(size(speciesCell));
    for k = 1:numel(speciesCell)
        nSnap = nCell{k};
        nSpeciesSnap = speciesCell{k};
        phiSnap = zeros(size(nSpeciesSnap));
        active = nSnap > 0;
        for s = 1:size(nSpeciesSnap,3)
            component = nSpeciesSnap(:,:,s);
            tmp = zeros(size(nSnap));
            tmp(active) = component(active)./nSnap(active);
            phiSnap(:,:,s) = tmp;
        end
        phiCell{k} = phiSnap;
    end
else
    error('snaps.phi or snaps.n_species is required for the species evolution figure.');
end

T = get_final_time(params, tSnap);
nt = numel(targetTimes);
snapData.t = zeros(1, nt);
snapData.tTarget = targetTimes(:)';
snapData.n = cell(1, nt);
snapData.phi = cell(1, nt);

for q = 1:nt
    tq = targetTimes(q);
    useFinal = (q == nt) || (~isnan(T) && abs(tq - T) <= 1e-12);

    if useFinal
        nSnap = nFinal;
        phiSnap = phiFinal;
        if isnan(T)
            tUsed = tq;
        else
            tUsed = T;
            if abs(T - tq) > 1e-8
                warning('Requested final panel t=%.6g, using final data at T=%.6g.', ...
                    tq, T);
            end
        end
    else
        [~, idx] = min(abs(tSnap - tq));
        idx = min(idx, min(numel(nCell), numel(phiCell)));
        nSnap = nCell{idx};
        phiSnap = phiCell{idx};
        tUsed = tSnap(idx);
        if abs(tUsed - tq) > 1e-8
            warning('Requested snapshot t=%.6g, using nearest saved snapshot t=%.6g.', ...
                tq, tUsed);
        end
    end

    snapData.t(q) = tUsed;
    snapData.n{q} = nSnap;
    snapData.phi{q} = phiSnap;
end
end

function cellData = snapshot_field_to_cell(A, spatialSize, fieldName)
if iscell(A)
    cellData = A(:);
    return;
end

if ndims(A) == 3 && size(A,1) == spatialSize(1) && size(A,2) == spatialSize(2)
    nt = size(A,3);
    cellData = cell(nt,1);
    for k = 1:nt
        cellData{k} = A(:,:,k);
    end
elseif ndims(A) == 3 && size(A,2) == spatialSize(1) && size(A,3) == spatialSize(2)
    nt = size(A,1);
    cellData = cell(nt,1);
    for k = 1:nt
        cellData{k} = squeeze(A(k,:,:));
    end
elseif ndims(A) == 4 && size(A,1) == spatialSize(1) && size(A,2) == spatialSize(2)
    nt = size(A,4);
    cellData = cell(nt,1);
    for k = 1:nt
        cellData{k} = A(:,:,:,k);
    end
elseif ndims(A) == 4 && size(A,2) == spatialSize(1) && size(A,3) == spatialSize(2)
    nt = size(A,1);
    cellData = cell(nt,1);
    for k = 1:nt
        cellData{k} = squeeze(A(k,:,:,:));
    end
else
    error('Unrecognized %s format.', fieldName);
end
end

function format_axes(ax, plotOpt)
set(ax, 'FontName', plotOpt.fontName, ...
    'FontSize', plotOpt.fontSizeAxis, ...
    'LineWidth', plotOpt.axisLineWidth, ...
    'Box', plotOpt.boxOn, ...
    'TickDir', plotOpt.tickDir, ...
    'Layer', 'top');
try
    set(ax, 'TickLabelInterpreter', plotOpt.interpreter);
catch
end
end

function apply_grid(ax, plotOpt)
if plotOpt.gridOn
    grid(ax, 'on');
else
    grid(ax, 'off');
end
end

function T = get_final_time(params, tVec)

T = NaN;

if isstruct(params)
    if isfield(params, 'T')
        T = params.T;
    elseif isfield(params, 'tFinal')
        T = params.tFinal;
    elseif isfield(params, 'tfinal')
        T = params.tfinal;
    elseif isfield(params, 'finalTime')
        T = params.finalTime;
    end
end

if isnan(T) && nargin >= 2 && ~isempty(tVec)
    T = tVec(end);
end

end

%% ========================================================================
function growthCoeff = get_growth_coeff(params, m)
% Return coefficients a_i in H_i(p)=a_i(1-p).
%
% For Example 4, the default is [8,4,1.5].  If your simulation stores these
% coefficients under a different field name, add that field below.

growthCoeff = [];

if isstruct(params)
    if isfield(params, 'growthCoeff')
        growthCoeff = params.growthCoeff;
    elseif isfield(params, 'growth_coeff')
        growthCoeff = params.growth_coeff;
    elseif isfield(params, 'aCoeff')
        growthCoeff = params.aCoeff;
    elseif isfield(params, 'a')
        growthCoeff = params.a;
    end
end

if isempty(growthCoeff)
    if m == 3
        growthCoeff = [8, 4, 1.5];
    else
        growthCoeff = ones(1,m);
    end
end

growthCoeff = growthCoeff(:)';

if numel(growthCoeff) < m
    growthCoeff = [growthCoeff, repmat(growthCoeff(end), 1, m-numel(growthCoeff))];
elseif numel(growthCoeff) > m
    growthCoeff = growthCoeff(1:m);
end

end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function save_all_formats(fig, baseName, savePdf, savePng, saveEps, pngResolution)
% Save selected formats. exportgraphics is used when available.

set(fig, 'Color', 'w');
set(fig, 'PaperPositionMode', 'auto');

if savePdf
    try
        exportgraphics(fig, [baseName '.pdf'], 'ContentType', 'vector');
    catch
        print(fig, [baseName '.pdf'], '-dpdf', '-painters');
    end
end

if savePng
    try
        exportgraphics(fig, [baseName '.png'], 'Resolution', pngResolution);
    catch
        print(fig, [baseName '.png'], '-dpng', sprintf('-r%d', pngResolution));
    end
end

if saveEps
    try
        print(fig, [baseName '.eps'], '-depsc2', '-painters');
    catch
        warning('Could not save EPS file: %s.eps', baseName);
    end
end

end
