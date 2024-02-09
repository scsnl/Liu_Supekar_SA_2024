function plot_explained_ratio(PCTVAR,dim,output_name)

plot(1:dim,100*PCTVAR(2,1:dim),'-o','LineWidth',2,'Color',[140/255,0,0]);
set(gca,'Fontsize',20);
set(gca,'LineWidth',2);
xlabel('Number of PLS components','FontSize',20);
ylabel('Explained ratio (%)','FontSize',20);
yticks([0:10:60])
grid off
box off

set(gcf,'position',[100,100,300,250])
print(gcf,'-dtiff',output_name);
close all

end