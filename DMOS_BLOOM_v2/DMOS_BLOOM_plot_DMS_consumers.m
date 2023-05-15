function [] = DMOS_BLOOM_plot_DMS_consumers(day_Meth,Meth_free,Meth_part,ndays,color,faceColor,lineType,fignum)
%DMOS_BLOOM_plot_DMS_consumers
%   Load relative abundances of Methylophaga from input file

figure(fignum)
subplot(2,2,1)
hold on
for b = 1:8
    plot(day_Meth(b,:),Meth_free(b,:),lineType{b},'MarkerSize',12,'Color',color(b,:),'MarkerFaceColor',faceColor(b,:),'LineWidth',2);
end
xlim([0 ndays])
ylim([0 0.016])
ymax = 0.016;
ymin = 0;
plot([ 8  8],[ymin ymax],'k--','LineWidth',5)
plot([13 13],[ymin ymax],'k--','LineWidth',5)
plot([18 18],[ymin ymax],'k--','LineWidth',5)
plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',12,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',12,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',8)
ylabel('Observed Methylophaga relative abundance','fontsize',8)
grid on
subplot(2,2,2)
hold on
for b = 1:8
    plot(day_Meth(b,:),Meth_part(b,:),lineType{b},'MarkerSize',12,'Color',color(b,:),'MarkerFaceColor',faceColor(b,:),'LineWidth',2);
end
xlim([0 ndays])
ylim([0 0.06])
ymax = 0.06;
ymin = 0;
plot([ 8  8],[ymin ymax],'k--','LineWidth',5)
plot([13 13],[ymin ymax],'k--','LineWidth',5)
plot([18 18],[ymin ymax],'k--','LineWidth',5)
plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',12,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',12,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',8)
ylabel('Observed Methylophaga relative abundance','fontsize',8)
grid on

return