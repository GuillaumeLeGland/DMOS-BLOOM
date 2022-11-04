function [] = DMOS_BLOOM_plot_panel(data,day,color,faceColor,lineType,ndays,ylab)

hold on

num_bag = size(data,1);

for b = 1:num_bag
    id = ~isnan(data(b,:));
    plot(day(b,id),data(b,id),lineType{b},'MarkerSize',5,'Color',color(b,:),'MarkerFaceColor',faceColor(b,:),'LineWidth',2);
end

yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',2)
plot([13 13],[0 ymax],'k--','LineWidth',2)
plot([18 18],[0 ymax],'k--','LineWidth',2)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])

xlim([0 ndays])
xlabel('Time [days]')
ylabel(ylab)

grid on

end