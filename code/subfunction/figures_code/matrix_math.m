function matrix_math(x,Ylabel,pic_size,output_path)


x = x([(size(x,1):-1:1)],:);
x = [x zeros(size(x,1),1)];
x = [x;zeros(1,size(x,2))];

figure;

m=pcolor(x);
m.LineWidth = 3;

J = customcolormap(linspace(0,1,2), {'#000000','#ffffff'});
h=colormap(J);

l=colorbar('southoutside','FontSize',15,'Box','on','LineWidth',1.5,'Ticks',[0:0.5:1.5],'TickLabels',[0:0.5:1.5]);

xlabel('Brain measure weight (five-percentile bins)');
yticks(1.5:1:(size(x,1)-1+0.5));

Ylabel = Ylabel([(size(x,1)-1):-1:1]);
yticklabels(Ylabel);
xticks(0);
set(gca,'FontSize',16);
caxis([0 1.7])

set(gcf,'position',[50,50,pic_size]);
print(gcf,'-dtiff',output_path);
close all
end