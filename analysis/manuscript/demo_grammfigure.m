
%% demo_grammfigure

x = sessiondata.treatment;
y = sessiondata.peak_force;
g=gramm('x',x,'y',y,'color',x);
% g.stat_violin('fill','transparent');
g.set_order_options('x',{'clonidine','placebo','atomoxetine'},...
                    'color',{'clonidine','placebo','atomoxetine'});
g.set_color_options('map',[0.54 0.27 0.07;...
                        [1 1 1]*0.5;...
                        1 0.5 0]);
g.geom_jitter('width',0.4);
g.set_color_options('map',[0.54 0.27 0.07;...
                        [1 1 1]*0.5;...
                        1 0.5 0],...
                    'chroma',20);
g.stat_boxplot();
% g.stat_summary('geom','bar','type','sem');
% g.stat_summary('geom','black_errorbar','type','sem');

g.draw();
for i=1:numel(g.results.stat_boxplot)
   g.results.stat_boxplot(i).box_handle.FaceAlpha = 0.75;
end