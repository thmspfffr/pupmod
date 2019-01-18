function f = figure_white(xSize,ySize)

f=figure;

set(gcf,'color','w');

% xSize = 20; ySize = 20;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')


