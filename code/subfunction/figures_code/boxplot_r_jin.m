function boxplot_r_jin(Perm_r,corr_r,Ylabel,pic_size,output_path)

figure;
[B,I] = sort(corr_r(:,1),'ascend');
Perm_r_new = Perm_r(I,:);
Ylabel_new = Ylabel(I)

boxplot(Perm_r_new','Orientation','horizontal','Whisker',3,'color',[137 137 137]/255)
set(gca,'FontSize',18);
set(gca,'LineWidth',3);
set(findobj(gca,'type','line'),'linew',3)
box off
hold on
xlim([-0.03 0.4])
for i=1:length(B)
 plot(B(i),i, 'o', 'MarkerSize', 15, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'yellow','LineWidth',2);
end

yticklabels(Ylabel_new);
set(gcf,'position',[50,50,pic_size]);
print(gcf,'-dtiff',output_path);
close all
end