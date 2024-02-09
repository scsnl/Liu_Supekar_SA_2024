function histogram_Jin(data, color, x_label, y_label, pic_title,output_path)

% h=histfit(data);
% h=histogram(data,'LineWidth',2,'FaceColor',color./255,'BinWidth',5);
 h=histogram(data,'LineWidth',2,'FaceColor',color./255);
xlabel(x_label);
ylabel(y_label);
title(pic_title);
set(gca,'FontSize',20,'LineWidth',3);
ymax_value=round(max(h.Values)+max(h.Values)./8);
ylim([0 ymax_value]);
% xmax_value=round(max(h.BinEdges)+range(h.BinEdges)./10);
% xmin_value=round(min(h.BinEdges)-range(h.BinEdges)./10);
% xlim([xmin_value xmax_value]);
% xlim([-10 60]);

set(gcf,'position',[250,250,370,280]);
print(gcf,'-dtiff','-r300',output_path);
% close

end