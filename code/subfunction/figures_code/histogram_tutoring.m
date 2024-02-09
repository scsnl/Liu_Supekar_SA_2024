function histogram_tutoring(data, color, x_label, y_label, pic_title,output_path)

% h=histfit(data);
h=histogram(data,5,'LineWidth',2,'FaceColor',color./255,'EdgeColor',[1 1 1]);
xlabel(x_label);
ylabel(y_label);
title(pic_title);
set(gca,'FontSize',20,'LineWidth',2);
box off
ymax_value=round(max(h.Values)+max(h.Values)./8);
ylim([0 ymax_value]);
xmax_value=round(max(h.BinEdges)+range(h.BinEdges)./10);
xmin_value=round(min(h.BinEdges)-range(h.BinEdges)./10);
xlim([xmin_value xmax_value]);

set(gcf,'position',[250,250,280,150]);
print(gcf,'-dtiff','-r300',output_path);
close

end