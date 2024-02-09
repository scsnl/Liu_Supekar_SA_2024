function beh_corr_scatter_NKI(x,y,Xlabel1,Xlabel2,Xlabel3,Ylabel,pic_title,output_path)


dotcolor = [231 223 202]./255;
linecolor = [182/255 180/255 194/255];
shadecolor = [154/255 154/255 154/255];

[xData, yData] = prepareCurveData(x, y);
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );
hold on

xFit = linspace(min(xData),max(xData),100);
yPredict = predint(fitresult,xFit,0.95,'functional','off');
fy = cat(2,yPredict(:,2)',flip(yPredict(:,1),1)')';
fx  = cat(2,xFit,flip(xFit',1)')';
fill(fx,fy,shadecolor,'EdgeAlpha',0,'FaceAlpha',0.25);
h1=plot( fitresult, xData, yData);
set(h1(1),'Marker','o','MarkerSize',10,'MarkerFaceColor',dotcolor,'MarkerEdgeColor','k')
set(h1(2),'LineWidth',4,'Color',shadecolor)

box off
set(gca,'FontSize',20);
set(gca,'LineWidth',3);
d1=(max([x])-min([x]))/6;
set(gca,'XLim',[min([x])-d1 max([x])+d1]);
d2=(max([y])-min([y]))/6;
set(gca,'YLim',[min([y])-d2 max([y])+d2]);

hold off
legend off

ylabel(Ylabel)
xlabel({Xlabel1;Xlabel2;Xlabel3})
title(pic_title)
% legend('Fitted','ASD','TD','Location','NorthEastOutside')  
set(gcf,'position',[50,50,350,390])
print(gcf,'-dtiff',output_path);
% close all
end