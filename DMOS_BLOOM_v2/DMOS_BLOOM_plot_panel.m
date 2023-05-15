function [] = DMOS_BLOOM_plot_panel(data,day,color,faceColor,lineType,ndays,ylab)

hold on

num_bag = size(data,1);

order = [1,2,3,5,6,4,7,8]; % Plot bag 4 after bags 5 and 6

for nb = 1:num_bag
    b = order(nb);
    id = ~isnan(data(b,:));
    plot(day(b,id),data(b,id),lineType{b},'MarkerSize',12,'Color',color(b,:),'MarkerFaceColor',faceColor(b,:),'LineWidth',2);
end

yl = ylim;
%ymax = yl(2);
odm = floor(log10(0.5*max(data(:)))); % Order of magnitude (2, 20, 200 ...) 
ymax = (10^odm)*ceil(max((10^(-odm))*data(:)));
plot([ 8  8],[0 ymax],'k--','LineWidth',5)
plot([13 13],[0 ymax],'k--','LineWidth',5)
plot([18 18],[0 ymax],'k--','LineWidth',5)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',12,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',12,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])

xlim([0 ndays])
xlabel('Time [days]')
ylabel(ylab)

grid on

end