function line_pre_post(pre,post,doc_color,mean_color,Ylim_range,YTick,YTicklabel, XTicklabel, pic_title,y_label, output_path)

data_mean = [mean(pre) mean(post)];
data = [pre post];

% pics for dis
c=plot(data','LineWidth',3,'Color',doc_color) ;
hold on
d=plot(data_mean','LineWidth',4,'Color',mean_color) ;
set(gca,'YLim',Ylim_range);
set(gca, 'YTick', YTick);
set(gca, 'YTickLabel', YTicklabel);
set(gca,'FontSize',18);
set(gca,'LineWidth',2);
set(gca,'XLim',[0.8 2.2]);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel' ,XTicklabel);
title(pic_title);
ylabel(y_label);
box off
set(gcf,'position',[100,100,280,300]) ;
print(gcf,'-dtiff',output_path);
close all
end