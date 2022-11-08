function [] = plotandsave_pip(pipmat, outputfile, muvals)

logmuvals = log10(muvals);
lenmu = length(muvals);
incr = (lenmu-1)/(logmuvals(end)-logmuvals(1));
ticks = [[1] [1 2] * incr lenmu];

contourf(pipmat-1,[-inf 0 inf]);
colormap([0 0 0; 1 1 1]);
set(gca, 'Units', 'inches', 'Position', [1 1 2 2]); axis image;
set(gca,'xtick', ticks,'ytick',ticks,'xticklabel',[],'yticklabel',[])
set(gcf, 'Units', 'inches', 'Position', [0 0 4 4]);
%set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.0 2.0]); 

exportgraphics(gcf, outputfile, 'BackgroundColor', 'none', 'ContentType', 'vector')