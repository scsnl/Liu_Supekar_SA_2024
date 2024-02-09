function beh_corr_scatter(x,y,Xlabel1,Xlabel2,Xlabel3,Ylabel,dotcolor,linecolor,pic_title,output_path)


% dotcolor = [170/255 201/255 206/255];
% linecolor = [182/255 180/255 194/255];
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
fill(fx,fy,shadecolor,'EdgeAlpha',0,'FaceAlpha',0.3);
h1=plot( fitresult, xData, yData);
set(h1(1),'Marker','.','MarkerSize',25,'Color',dotcolor)
set(h1(2),'LineWidth',4,'Color',linecolor)

% [p] = polyfit(x,y,1);
%   f = polyval(p,x); 
%  plot(x,y,'.','MarkerEdgeColor',[157/255,211/255,250/255],'MarkerSize',35)
%  hold on
% 
% plot(x,f,'k-','LineWidth',4,'MarkerSize',20) 
box off
set(gca,'FontSize',18);
set(gca,'LineWidth',2);
% plot(x(group==1),y(group==1),'.','MarkerEdgeColor',[250/255,128/255,114/255],'MarkerSize',30)
% plot(x(group==2),y(group==2),'.','MarkerEdgeColor',[61/255,89/255,171/255],'MarkerSize',30)
% set(gca,'XLim',[0.5 5.5]);
% set(gca, 'XTick', [1:5])
% set(gca, 'XTickLabel' ,{'Day1','Day2','Day3','Day4','Day5'})
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
set(gcf,'position',[50,50,300,300])
print(gcf,'-dtiff',output_path);
close all
end