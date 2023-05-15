function [] = DMOS_BLOOM_plot_model(zoo_mod,Chl_mod,pico_mod,nano_mod,nano1_mod,nano2_mod,Ehux_mod,bac_mod,...
                                    DIN_mod,ON_mod,DMSPp_mod,DMSPd_mod,DMS_mod,FDC_mod,...
                                    obsPlankton,cil_obs,Chl_obs,pico_obs,nano_obs,Ehux_obs,bac_obs,NO3_obs,...
                                    obsSulfur,DMSPp_obs,DMSPd_obs,DMS_obs,PON_obs,offset,tspan,day,color,fignum)


%===============================================================================
%...............................................................................
%===============================================================================

% For figure clarity, the model concentrations of nanophytoplankton groups
% 1 and 2 (nano1_mod and nano2_mod) are not shown in the figures.
% This can be changed by remmoving comments in subplot 2 of Fig. 2.

% Color for bags 1-2-3-5-6: [0 0.8 0.6]
figure(fignum(1))
subplot(2,2,1);
hold on
plot(tspan,zoo_mod(:,1),'k')
plot(tspan,zoo_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,zoo_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,zoo_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
plot(day(1,:),obsPlankton(1,:,4)','k.','MarkerSize',15);
plot(day(1,:),obsPlankton(2,:,4)','r.','MarkerSize',15);
plot(day(1,:),cil_obs(7,:),'.','MarkerSize',15,'Color',color(7,:));
plot(day(1,:),cil_obs(4,:),'.','MarkerSize',15,'Color',color(4,:));
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',15)
ylabel('Zooplankton concentration [mmolN m^{-3}]','fontsize',15)
grid on
subplot(2,2,2)
hold on
plot(tspan,Chl_mod(:,1),'k')
plot(tspan,Chl_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,Chl_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,Chl_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
for b = 1:8
    plot(day(b,:),Chl_obs(b,:),'.','MarkerSize',15,'Color',color(b,:));
end
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',12)
ylabel('Model chlorophyll concentration [mgChl m^{-3}]','fontsize',12)
grid on

%===============================================================================
%...............................................................................
%===============================================================================

figure(fignum(2))
subplot(2,3,1);
hold on
plot(day(1,:),obsPlankton(1,:,1)','kd','MarkerSize',9,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsPlankton(2,:,1)','p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),pico_obs(7,:),'s','MarkerSize',9,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),pico_obs(4,:),'o','MarkerSize',9,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,pico_mod(:,1),'k','LineWidth',3)
plot(tspan,pico_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,pico_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,pico_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
%errorbar(day(1,:),obs_plankton(2,:,1)',obs_plankton_tol(2,:,1)',obs_plankton_tol(2,:,1)',opt{:})
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',15)
ylabel('Picophytoplankton concentration [mmolN m^{-3}]','fontsize',15)
grid on
subplot(2,3,2);
hold on
plot(day(1,:),obsPlankton(1,:,2)','kd','MarkerSize',9,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsPlankton(2,:,2)','p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),nano_obs(7,:),'s','MarkerSize',9,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),nano_obs(4,:),'o','MarkerSize',9,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,nano_mod(:,1),'k','LineWidth',3)
%plot(tspan,nano1_mod(:,1),'k--','LineWidth',3)
%plot(tspan,nano2_mod(:,1),'k--','LineWidth',3)
plot(tspan,nano_mod(:,4),'Color',color(7,:),'LineWidth',3)
%plot(tspan,nano1_mod(:,4),'--','Color',color(7,:),'LineWidth',3)
%plot(tspan,nano2_mod(:,4),'--','Color',color(7,:),'LineWidth',3)
plot(tspan,nano_mod(:,3),'Color',color(4,:),'LineWidth',3)
%plot(tspan,nano1_mod(:,3),'--','Color',color(4,:),'LineWidth',3)
%plot(tspan,nano2_mod(:,3),'--','Color',color(4,:),'LineWidth',3)
plot(tspan,nano_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
%plot(tspan,nano1_mod(:,2),'--','Color',[0.5 0.5 0.5],'LineWidth',3)
%plot(tspan,nano2_mod(:,2),'--','Color',[0.5 0.5 0.5],'LineWidth',3)
%errorbar(day(1,:),obs_plankton(2,:,2)',obs_plankton_tol(2,:,2)',obs_plankton_tol(2,:,2)',opt{:});
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',15)
ylabel('nanophytoplankton concentration [mmolN m^{-3}]','fontsize',15)
grid on
subplot(2,3,3);
hold on
plot(day(1,:),obsPlankton(1,:,3)','kd','MarkerSize',9,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsPlankton(2,:,3)','p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),Ehux_obs(7,:),'s','MarkerSize',9,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),Ehux_obs(4,:),'o','MarkerSize',9,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,Ehux_mod(:,1),'k','LineWidth',3)
plot(tspan,Ehux_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,Ehux_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,Ehux_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
%errorbar(day(1,:),obs_plankton(2,:,3)',obs_plankton_tol(2,:,3)',obs_plankton_tol(2,:,3)',opt{:})
%yl = ylim;
ymax = 14.0; %ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',15)
ylabel('Emiliania huxleyi concentration [mmolN m^{-3}]','fontsize',15)
grid on
subplot(2,3,4);
hold on
plot(day(1,:),obsPlankton(1,:,5)','kd','MarkerSize',9,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsPlankton(2,:,5)','p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),bac_obs(7,:),'s','MarkerSize',9,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),bac_obs(4,:),'o','MarkerSize',9,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,bac_mod(:,1),'k','LineWidth',3)
plot(tspan,bac_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,bac_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,bac_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
%errorbar(day(1,:),obs_plankton(2,:,5)',obs_plankton_tol(2,:,5)',obs_plankton_tol(2,:,5)',opt{:})
%yl = ylim;
ymax = 1.6; %ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',15)
ylabel('Bacteria concentration [mmolN m^{-3}]','fontsize',15)
grid on
subplot(2,3,5);
hold on
plot(day(1,:),obsPlankton(1,:,6)' - offset,'kd','MarkerSize',9,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsPlankton(2,:,6)' - offset,'p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),NO3_obs(7,:),'s','MarkerSize',9,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),NO3_obs(4,:),'o','MarkerSize',9,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,DIN_mod(:,1),'k','LineWidth',3)
plot(tspan,DIN_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,DIN_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,DIN_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',15)
ylabel('DIN concentration [mmolN m^{-3}]','fontsize',15)
grid on
subplot(2,3,6);
hold on
plot(day(1,:),obsPlankton(1,:,7)' - offset,'kd','MarkerSize',9,'MarkerFaceColor',[1 1 1],'LineWidth',2);
plot(day(1,:),obsPlankton(2,:,7)' - offset,'p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[1 1 1],'Markersize',15,'LineWidth',2);
plot(day(1,:),PON_obs(7,:),'s','MarkerSize',9,'Color',color(7,:),'MarkerFaceColor',[1 1 1],'LineWidth',2);
plot(day(1,:),PON_obs(4,:),'o','MarkerSize',9,'Color',color(4,:),'MarkerFaceColor',[1 1 1],'LineWidth',2);
plot(tspan,ON_mod(:,1),'k','LineWidth',3)
plot(tspan,ON_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,ON_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,ON_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
ylim([0 14])
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',9,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
xlabel('Time [days]','fontsize',15)
ylabel('ON concentration [mmolN m^{-3}]','fontsize',15)
grid on

%===============================================================================
%...............................................................................
%===============================================================================

figure(fignum(3))
subplot(4,4,1);
hold on
ylim([0 800])
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(day(1,:),obsSulfur(1,:,3)','kd','MarkerSize',6,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsSulfur(2,:,3)','p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),DMSPp_obs(7,:),'s','MarkerSize',6,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),DMSPp_obs(4,:),'o','MarkerSize',6,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,DMSPp_mod(:,1),'k','LineWidth',3)
plot(tspan,DMSPp_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,DMSPp_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,DMSPp_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('Time [days]','fontsize',8)
ylabel('Particulate DMSP concentration [umolS m^{-3}]','fontsize',8)
grid on
subplot(4,4,5);
hold on
ylim([0 100])
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(day(1,:),obsSulfur(1,:,2)','kd','MarkerSize',6,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsSulfur(2,:,2)','p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),DMSPd_obs(7,:),'s','MarkerSize',6,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),DMSPd_obs(4,:),'o','MarkerSize',6,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,DMSPd_mod(:,1),'k','LineWidth',3)
plot(tspan,DMSPd_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,DMSPd_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,DMSPd_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('Time [days]','fontsize',8)
ylabel('Dissolved DMSP concentration [umolS m^{-3}]','fontsize',8)
grid on
subplot(4,4,9);
hold on
ylim([0 140])
set(gca,'ytick',0:20:140)
yl = ylim;
ymax = yl(2);
plot([ 8  8],[0 ymax],'k--','LineWidth',3)
plot([13 13],[0 ymax],'k--','LineWidth',3)
plot([18 18],[0 ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,0.975*ymax*ones(1,8),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,0.975*ymax*ones(1,5),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(day(1,:),obsSulfur(1,:,1)','kd','MarkerSize',6,'MarkerFaceColor',[0 0 0],'LineWidth',2);
plot(day(1,:),obsSulfur(2,:,1)','p','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Markersize',15);
plot(day(1,:),DMS_obs(7,:),'s','MarkerSize',6,'Color',color(7,:),'MarkerFaceColor',color(7,:),'LineWidth',2);
plot(day(1,:),DMS_obs(4,:),'o','MarkerSize',6,'Color',color(4,:),'MarkerFaceColor',color(4,:),'LineWidth',2);
plot(tspan,DMS_mod(:,1),'k','LineWidth',3)
plot(tspan,DMS_mod(:,4),'Color',color(7,:),'LineWidth',3)
plot(tspan,DMS_mod(:,3),'Color',color(4,:),'LineWidth',3)
plot(tspan,DMS_mod(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('Time [days]','fontsize',8)
ylabel('DMS concentration [umolS m^{-3}]','fontsize',7)
grid on
subplot(4,4,13);
hold on
ylim([0 0.04])
yl = ylim;
ymax = yl(2);
ymin = yl(1);
plot([ 8  8],[ymin ymax],'k--','LineWidth',3)
plot([13 13],[ymin ymax],'k--','LineWidth',3)
plot([18 18],[ymin ymax],'k--','LineWidth',3)
plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',6,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
plot(tspan,exp(FDC_mod(:,1))./(1+exp(FDC_mod(:,1))),'k','LineWidth',3)
plot(tspan,exp(FDC_mod(:,4))./(1+exp(FDC_mod(:,4))),'Color',color(7,:),'LineWidth',3)
plot(tspan,exp(FDC_mod(:,3))./(1+exp(FDC_mod(:,3))),'Color',color(4,:),'LineWidth',3)
plot(tspan,exp(FDC_mod(:,2))./(1+exp(FDC_mod(:,2))),'Color',[0.5 0.5 0.5],'LineWidth',3)
xlabel('Time [days]','fontsize',8)
ylabel('Fraction of DMS-consuming bacteria [no unit]','fontsize',8)
grid on

%===============================================================================
%...............................................................................
%===============================================================================


end