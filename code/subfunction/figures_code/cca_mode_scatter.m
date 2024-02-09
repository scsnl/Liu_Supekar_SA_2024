function cca_mode_scatter(x,y,Xlabel1,Ylabel1,output_path)

figure('position',[100,100,430,420],'PaperPositionMode','auto')
scatter(x, y, 120, 'MarkerEdgeColor','k','MarkerFaceColor',[217,237,232]./255,'LineWidth',1.5);
hold on
set(gca,'FontSize',25);
set(gca,'LineWidth',4);
xlabel(Xlabel1);
ylabel(Ylabel1);
xticks([-3 -2 -1 0 1 2 3]);
yticks([-3 -2 -1 0 1 2 3]);

d1=(max([x])-min([x]))/10;
set(gca,'XLim',[min([x])-d1 max([x])+d1]);
d2=(max([y])-min([y]))/10;
set(gca,'YLim',[min([y])-d2 max([y])+d2]);
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[233 233 233]./255);
set(0,'DefaultFigureVisible', 'on') 
print(gcf,'-dtiffn','-r300',output_path);
% close all

end