function bar_GOenrich(data,Y_label,title_label,top_num,output_name)

figure;
% J = customcolormap(linspace(0,1,10), {'#cf1619','#d87323','#e8b52c','#f9fc35','#4fa019','#85fb27','#6a70b1','#422450','#430e32','#000000'});
% colormap(J)
J = colormap;
[ranked_b ind]= sort(data.b,'ascend');
h=barh(ranked_b(1:top_num),'LineWidth',2);
box off
set(gca,'Fontsize',20)
set(gca,'LineWidth',2)
xlabel('Number of overlapped genes','FontSize',20);
ylabel(Y_label,'FontSize',20);
ranked_go_label = data.Description(ind);
ranked_q=data.FDRQ_value(ind);
log_q=-log10(ranked_q);
a = 1;
b = 64;
Ymax = max(log_q);
Ymin = min(log_q);
k = (b-a)/(Ymax-Ymin);
scale_q=a+k*(log_q-Ymin);
yticks([]);
title(title_label,'FontSize',20);

h.FaceColor='flat';
for i=1:top_num
    color_ind=round(scale_q(i));
    if color_ind <=1
        color_ind = 1;
    elseif color_ind >64
        color_ind = 64;
    end
  h.CData(i,:) = J(color_ind,:); 
   text(max(ranked_b(1:top_num))*0.02,i,['GO:' ranked_go_label{i}],'FontSize',15);
end

c=colorbar('southoutside','FontSize',20,'Box','off','Ticks',[0:0.2:1],'TickLabels',[4:1:9]);
c.Label.String = '-log10(q)';
c.Label.FontSize = 15;
set(gcf,'position',[100,100,555,500]);
print(gcf,'-dtiff',output_name);
 close all

end
