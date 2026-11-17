function import_manuscript_curves(manuscript_dir,data_root)
%IMPORT_MANUSCRIPT_CURVES Extract numeric plotting arrays from source FIG files.
% This import is optional provenance tooling, not a simulation or redraw step.
% The saved MAT files contain numeric arrays and text only, never graphics
% objects, callbacks, function handles, source code, or machine paths.
figures={{'temporal_convergence.fig','spatial_convergence.fig'}, ...
    {'total_structure.fig'},{'two_three_species3.fig'}};
for k=1:3
    dest=fullfile(data_root,sprintf('Figure%d',k));
    if ~exist(dest,'dir'),mkdir(dest);end
    plot_data=struct('figure_number',k,'panels',[]);
    for j=1:numel(figures{k})
        f=openfig(fullfile(manuscript_dir,figures{k}{j}),'invisible');
        cleaner=onCleanup(@()close(f));
        axs=findall(f,'Type','axes');[~,order]=sort(arrayfun(@(a)a.Position(1),axs));axs=axs(order);
        for a=reshape(axs,1,[])
            p=extract_panel(a,f);
            if isempty(plot_data.panels),plot_data.panels=p;else,plot_data.panels(end+1)=p;end
        end
        clear cleaner
    end
    save(fullfile(dest,'plot_data.mat'),'plot_data','-v7');
    for p=1:numel(plot_data.panels)
        for c=1:numel(plot_data.panels(p).curves)
            z=plot_data.panels(p).curves(c);
            if strcmp(z.role,'data')
                writematrix([z.x(:),z.y(:)],fullfile(dest,sprintf('panel%02d_curve%02d.csv',p,c)));
            end
        end
    end
end
end
function p=extract_panel(a,f)
p=struct('title',plain(a.Title.String),'xlabel',plain(a.XLabel.String), ...
    'ylabel_left','','ylabel_right','','xlim',a.XLim,'ylim_left',[], ...
    'ylim_right',[],'xscale',a.XScale,'yscale',a.YScale, ...
    'font_size',a.FontSize,'position',a.Position,'curves',[], ...
    'legend_indices',[],'legend_labels',{{}},'legend_location','northeast', ...
    'legend_font_size',a.FontSize,'legend_columns',1,'texts',[]);
if numel(a.YAxis)==2
    yyaxis(a,'left');p.ylabel_left=plain(a.YLabel.String);p.ylim_left=a.YLim;
    yyaxis(a,'right');p.ylabel_right=plain(a.YLabel.String);p.ylim_right=a.YLim;
else
    p.ylabel_left=plain(a.YLabel.String);p.ylim_left=a.YLim;
end
lines=findall(a,'Type','line');lines=flipud(lines(:));
for i=1:numel(lines)
    l=lines(i);side=1;
    if numel(a.YAxis)==2
        yyaxis(a,'right');if any(a.Children==l),side=2;end
    end
    name=plain(l.DisplayName);role='data';
    if isempty(name) && numel(l.XData)<=3
        if contains(p.title,'mass fractions'),role='legend';else,role='reference';end
    end
    marker_indices=[];if isprop(l,'MarkerIndices'),marker_indices=l.MarkerIndices;end
    c=struct('x',double(l.XData),'y',double(l.YData),'name',name, ...
        'side',side,'role',role,'color',l.Color,'line_style',l.LineStyle, ...
        'line_width',l.LineWidth,'marker',l.Marker,'marker_size',l.MarkerSize, ...
        'marker_indices',marker_indices);
    if isempty(p.curves),p.curves=c;else,p.curves(end+1)=c;end
end
legends=findall(f,'Type','legend');
for lg=reshape(legends,1,[])
    pc=lg.PlotChildren;indices=[];
    for h=reshape(pc,1,[])
        ix=find(lines==h,1);if ~isempty(ix),indices(end+1)=ix;end
    end
    if ~isempty(indices)
        p.legend_indices=indices;p.legend_labels=cellstr(string(lg.String));
        p.legend_location=lg.Location;if strcmp(p.legend_location,'none'),p.legend_location='northeast';end
        p.legend_font_size=lg.FontSize;p.legend_columns=lg.NumColumns;
    end
end
custom=a.Children(strcmp(get_types(a.Children),'text'));
for t=reshape(custom,1,[])
    if isempty(t.String),continue;end
    z=struct('text',plain(t.String),'position',t.Position,'units',t.Units, ...
        'color',t.Color,'font_size',t.FontSize,'interpreter',t.Interpreter, ...
        'horizontal_alignment',t.HorizontalAlignment, ...
        'vertical_alignment',t.VerticalAlignment,'font_weight',t.FontWeight);
    if isempty(p.texts),p.texts=z;else,p.texts(end+1)=z;end
end
end
function types=get_types(h)
types=cell(numel(h),1);for j=1:numel(h),types{j}=h(j).Type;end
end
function s=plain(s)
if iscell(s),s=strjoin(s,' ');end
s=char(string(s));
end
