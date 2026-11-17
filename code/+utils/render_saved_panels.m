function output_files=render_saved_panels(data_root,out_dir,number)
%RENDER_SAVED_PANELS Redraw paper curves from plain numerical MAT data.
% No source FIG file or numerical solver is used by this function.
if ~exist(out_dir,'dir'),mkdir(out_dir);end
v=load(fullfile(data_root,sprintf('Figure%d',number),'plot_data.mat'),'plot_data');
P=v.plot_data.panels;count=numel(P);
if count==2,sz=[1 1 13.4 5.8];else,sz=[1 1 18 5.25];end
f=figure('Visible','off','Color','w','Units','inches','Position',sz,'Renderer','painters');
cleanup=onCleanup(@()close(f));
tiledlayout(f,1,count,'Padding','compact','TileSpacing','compact');
for i=1:count
    p=P(i);ax=nexttile;hold(ax,'on');
    set(ax,'XScale',p.xscale,'YScale',p.yscale,'Box','on', ...
        'TickDir','in','LineWidth',1.0,'FontName','Times New Roman', ...
        'FontSize',p.font_size,'TickLabelInterpreter','latex');
    handles=gobjects(1,numel(p.curves));
    custom_mass_legend=contains(p.title,'mass fractions');
    for j=1:numel(p.curves)
        c=p.curves(j);
        if custom_mass_legend && strcmp(c.role,'legend'),continue;end
        if ~isempty(p.ylim_right)
            if c.side==2,yyaxis(ax,'right');else,yyaxis(ax,'left');end
        end
        handles(j)=plot(ax,c.x,c.y,'Color',c.color,'LineStyle',c.line_style, ...
            'LineWidth',c.line_width,'Marker',c.marker,'MarkerSize',c.marker_size, ...
            'MarkerFaceColor','w','DisplayName',c.name);
        if ~isempty(c.marker_indices),handles(j).MarkerIndices=c.marker_indices;end
        if ~strcmp(c.role,'data'),handles(j).HandleVisibility='off';end
    end
    xlim(ax,p.xlim);
    if ~isempty(p.ylim_right)
        yyaxis(ax,'left');ylim(ax,p.ylim_left);
        ylabel(ax,p.ylabel_left,'Interpreter','latex','FontSize',p.font_size+1);
        yyaxis(ax,'right');ylim(ax,p.ylim_right);
        ylabel(ax,p.ylabel_right,'Interpreter','latex','FontSize',p.font_size+1);
    else
        ylim(ax,p.ylim_left);ylabel(ax,p.ylabel_left,'Interpreter','latex','FontSize',p.font_size+1);
    end
    xlabel(ax,p.xlabel,'Interpreter','latex','FontSize',p.font_size+1);
    title(ax,p.title,'Interpreter','latex','FontSize',p.font_size,'FontWeight','normal');
    if ~isempty(p.legend_indices)
        legend(ax,handles(p.legend_indices),p.legend_labels,'Interpreter','latex', ...
            'Location',p.legend_location,'FontSize',p.legend_font_size, ...
            'NumColumns',p.legend_columns,'Box','on','Color','w');
    end
    if custom_mass_legend
        draw_mass_legend(ax,p);
    end
    for t=p.texts
        if custom_mass_legend,continue;end
        text(ax,t.position(1),t.position(2),t.text,'Units',t.units, ...
            'Color',t.color,'FontSize',t.font_size,'Interpreter',t.interpreter, ...
            'HorizontalAlignment',t.horizontal_alignment, ...
            'VerticalAlignment',t.vertical_alignment,'FontWeight',t.font_weight);
    end
end
stems={'Figure1','total_structure','two_three_species3'};stem=fullfile(out_dir,stems{number});
exportgraphics(f,[stem '.pdf'],'ContentType','vector');
exportgraphics(f,[stem '.png'],'Resolution',180);
output_files={[stem '.pdf'],[stem '.png']};
fprintf('Redrew Figure %d from saved numerical arrays: %s\n',number,stem);
end

function draw_mass_legend(ax,p)
% Display-only legend samples are separated from numerical curve data.
rectangle(ax,'Position',[.04 .69 1.92 .29],'FaceColor','w', ...
    'EdgeColor',[.5 .5 .5],'LineWidth',.75);
font=16;
text(ax,.07,.94,'two species','Units','normalized','FontName','Times New Roman','FontSize',font);
text(ax,.55,.94,'three species','Units','normalized','FontName','Times New Roman','FontSize',font);
colors=[0 .447 .741;.85 .325 .098;0 .55 1;1 .12 .12;.2 .62 .2];
for k=1:5
    ix=find(arrayfun(@(c)strcmp(c.role,'data') && norm(c.color-colors(k,:))<1e-6,p.curves),1);
    c=p.curves(ix);group=1+(k>=3);row=k;if group==2,row=k-2;end
    x0=.14;if group==2,x0=1.10;end
    y=.86-(row-1)*.07;
    h=plot(ax,[x0,x0+.14,x0+.28],[y,y,y], ...
        'Color',c.color,'LineStyle',c.line_style,'LineWidth',c.line_width, ...
        'Marker',c.marker,'MarkerSize',c.marker_size,'MarkerFaceColor','w', ...
        'HandleVisibility','off');
    if ~strcmp(c.marker,'none'),h.MarkerIndices=2;end
    text(ax,x0+.33,y,sprintf('$M_%d/M$',row),'Interpreter','latex', ...
        'FontSize',font,'VerticalAlignment','middle');
end
end
