function cca_mode_beh_weight(x,label,output_path)

figure('position',[100,100,350,450],'PaperPositionMode','auto')
% J = customcolormap(linspace(0,1,11), {'#ff8200','#fdad31','#e7b62b','#fcda70','#fefec2','#ffffff','#d0fefe','#b3f7fe','#7ce2fe','#33b0fc','#007bf3'});
% J = customcolormap(linspace(0,1,11), {'#7abace','#499ace','#3987cf','#465dbd','#43417f','#000000','#291e20','#49242d','#a52f44','#c93f3a','#f0894b'});
J = customcolormap(linspace(0,1,11), {'#f8b07d','#f0894b','#c93f3a','#a52f44','#49242d','#000000','#43417f','#465dbd','#3987cf','#499ace','#7abace'});

h=colormap(J);

box off
hold on
set(gca,'FontSize',18);
set(gca,'xtick',[0:0:0]);

if x < 0
    set(gca,'yLim',[-1.2,0]);
    set(gca,'ytick',[sort(x),0]);
    for i=1:length(x)
        text(0.05,x(i),[label{i}],'FontSize',20,'Color',h(size(h,1)-round(-x(i)*size(h,1))+1,:),'FontWeight','bold')
    end
end

if x > 0
    set(gca,'yLim',[0,1.2]);
    set(gca,'ytick',[0,sort(x)]);
    for i=1:length(x)
        text(0.05,x(i),[label{i}],'FontSize',20,'Color',h(round(x(i)*size(h,1))-1,:),'FontWeight','bold')
    end
end

set(gca,'LineWidth',3);
colorbar('southoutside','FontSize',18,'Box','off','LineWidth',1.5,'Ticks',[0:0.25:1],'TickLabels',[-1:0.5:1]);
print(gcf,'-dtiffn','-r300',output_path);
% close all

end