function DMOS_BLOOM_plot_fluxes(f,DMSPd_mod,DMS_mod,tspan,color,fignum)
% Plot DMOS_BLOOM fluxes

%========================================================================

nphy = 4;
nzoo = 5;

%--------------------------------------------------------------------------
% Plankton model fluxes
%--------------------------------------------------------------------------

figure(fignum(1))

    subplot(3,3,1)
    hold on
    for ji = 1:nphy
        plot(tspan,f.uptake_DIN_phy(ji,:,1),'k')
        plot(tspan,f.uptake_DIN_phy(ji,:,4),'Color',color(7,:))
        plot(tspan,f.uptake_DIN_phy(ji,:,3),'Color',color(4,:))
        plot(tspan,f.uptake_DIN_phy(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',10)
    ylabel('Phytoplankton DIN uptake [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on
    
    subplot(3,3,2)
    hold on
    for ji = 1:nzoo
        plot(tspan,f.graz_zoo(ji,:,1),'k')
        plot(tspan,f.graz_zoo(ji,:,4),'Color',color(7,:))
        plot(tspan,f.graz_zoo(ji,:,3),'Color',color(4,:))
        plot(tspan,f.graz_zoo(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',10)
    ylabel('Zooplankton grazing [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on
    
    subplot(3,3,3)
    hold on
    for ji = 1:nphy
        plot(tspan,f.mort_phy(ji,:,1),'k')
        plot(tspan,f.mort_phy(ji,:,4),'Color',color(7,:))
        plot(tspan,f.mort_phy(ji,:,3),'Color',color(4,:))
        plot(tspan,f.mort_phy(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',10)
    ylabel('Phytoplankton mortality fluxes [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on

    subplot(3,3,4)
    hold on
    plot(tspan,f.ON_to_bac(:,1),'k')
    plot(tspan,f.ON_to_bac(:,4),'Color',color(7,:))
    plot(tspan,f.ON_to_bac(:,3),'Color',color(4,:))
    plot(tspan,f.ON_to_bac(:,2),'Color',color(2,:))
    xlabel('Time [days]','fontsize',10)
    ylabel('Bacteria growth [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on
    
    subplot(3,3,5)
    hold on
    plot(tspan,f.remin_ON(:,1),'k')
    plot(tspan,f.remin_ON(:,4),'Color',color(7,:))
    plot(tspan,f.remin_ON(:,3),'Color',color(4,:))
    plot(tspan,f.remin_ON(:,2),'Color',color(2,:))
    xlabel('Time [days]','fontsize',10)
    ylabel('Bacterial ON remineralization [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on
    
    subplot(3,3,6)
    hold on
    plot(tspan,f.mort_bac(:,1),'k')
    plot(tspan,f.mort_bac(:,4),'Color',color(7,:))
    plot(tspan,f.mort_bac(:,3),'Color',color(4,:))
    plot(tspan,f.mort_bac(:,2),'Color',color(2,:))
    xlabel('Time [days]','fontsize',10)
    ylabel('Bacteria mortality fluxes [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on
    
    subplot(3,3,7)
    hold on
    for ji = 1:nzoo
        plot(tspan,f.prey_to_zoo(ji,:,1),'k')
        plot(tspan,f.prey_to_zoo(ji,:,4),'Color',color(7,:))
        plot(tspan,f.prey_to_zoo(ji,:,3),'Color',color(4,:))
        plot(tspan,f.prey_to_zoo(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',10)
    ylabel('Zooplankton growth [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on
    
    subplot(3,3,8)
    hold on
    for ji = 1:nzoo
        plot(tspan,f.exud(ji,:,1),'k')
        plot(tspan,f.exud(ji,:,4),'Color',color(7,:))
        plot(tspan,f.exud(ji,:,3),'Color',color(4,:))
        plot(tspan,f.exud(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',10)
    ylabel('Zooplankton exudation [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on
    
    subplot(3,3,9)
    hold on
    for ji = 1:nzoo
        plot(tspan,f.mort_zoo(ji,:,1),'k')
        plot(tspan,f.mort_zoo(ji,:,4),'Color',color(7,:))
        plot(tspan,f.mort_zoo(ji,:,3),'Color',color(4,:))
        plot(tspan,f.mort_zoo(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',10)
    ylabel('Zooplankton mortality [mmolN m^{-3} d^{-1}]','fontsize',10)
    grid on

%--------------------------------------------------------------------------
% Sulfur model fluxes
%--------------------------------------------------------------------------    

figure(fignum(2))

    subplot(4,4,1)
    hold on
    for ji = 1:nphy
        plot(tspan,f.leak_DMS(ji,:,1),'k')
        plot(tspan,f.leak_DMS(ji,:,4),'Color',color(7,:))
        plot(tspan,f.leak_DMS(ji,:,3),'Color',color(4,:))
        plot(tspan,f.leak_DMS(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',8)
    ylabel('Phytoplankton DMS leakage [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,2)
    hold on
    plot(tspan,f.bac_to_DMS(:,1),'k')
    plot(tspan,f.bac_to_DMS(:,4),'Color',color(7,:))
    plot(tspan,f.bac_to_DMS(:,3),'Color',color(4,:))
    plot(tspan,f.bac_to_DMS(:,2),'Color',color(2,:))
    xlabel('Time [days]','fontsize',8)
    ylabel('Bacterial DMS production [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,3)
    hold on
    for ji = 1:nphy
        plot(tspan,f.mort_to_DMS(ji,:,1),'k')
        plot(tspan,f.mort_to_DMS(ji,:,4),'Color',color(7,:))
        plot(tspan,f.mort_to_DMS(ji,:,3),'Color',color(4,:))
        plot(tspan,f.mort_to_DMS(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',8)
    ylabel('Mortality to DMS [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,4)
    hold on
    for ji = 1:nphy
        plot(tspan,f.graz_to_DMS(ji,:,1),'k')
        plot(tspan,f.graz_to_DMS(ji,:,4),'Color',color(7,:))
        plot(tspan,f.graz_to_DMS(ji,:,3),'Color',color(4,:))
        plot(tspan,f.graz_to_DMS(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',8)
    ylabel('Grazing to DMS [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,5)
    hold on
    for ji = 1:nphy
        plot(tspan,f.leak_DMSP(ji,:,1),'k')
        plot(tspan,f.leak_DMSP(ji,:,4),'Color',color(7,:))
        plot(tspan,f.leak_DMSP(ji,:,3),'Color',color(4,:))
        plot(tspan,f.leak_DMSP(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',8)
    ylabel('Phytoplankton DMSP leakage [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,7)
    hold on
    for ji = 1:nphy
        plot(tspan,f.mort_to_DMSP(ji,:,1),'k')
        plot(tspan,f.mort_to_DMSP(ji,:,4),'Color',color(7,:))
        plot(tspan,f.mort_to_DMSP(ji,:,3),'Color',color(4,:))
        plot(tspan,f.mort_to_DMSP(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',8)
    ylabel('Mortality to DMSP [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,8)
    hold on
    for ji = 1:nphy
        plot(tspan,f.graz_to_DMSP(ji,:,1),'k')
        plot(tspan,f.graz_to_DMSP(ji,:,4),'Color',color(7,:))
        plot(tspan,f.graz_to_DMSP(ji,:,3),'Color',color(4,:))
        plot(tspan,f.graz_to_DMSP(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',8)
    ylabel('Grazing to DMSP [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,9)
    hold on
    plot(tspan,f.Udms(:,1),'k')
    plot(tspan,f.Udms(:,4),'Color',color(7,:))
    plot(tspan,f.Udms(:,3),'Color',color(4,:))
    plot(tspan,f.Udms(:,2),'Color',color(2,:))
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS uptake by bacteria [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on

    subplot(4,4,12)
    hold on
    plot(tspan,f.DMS_sink(:,1),'k')
    plot(tspan,f.DMS_sink(:,4),'Color',color(7,:))
    plot(tspan,f.DMS_sink(:,3),'Color',color(4,:))
    plot(tspan,f.DMS_sink(:,2),'Color',color(2,:))
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sink [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
    subplot(4,4,13)
    hold on
    plot(tspan,f.Udmsp(:,1),'k')
    plot(tspan,f.Udmsp(:,4),'Color',color(7,:))
    plot(tspan,f.Udmsp(:,3),'Color',color(4,:))
    plot(tspan,f.Udmsp(:,2),'Color',color(2,:))
    xlabel('Time [days]','fontsize',8)
    ylabel('Bulk bacterial DMSP uptake [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on

    subplot(4,4,14)
    hold on
    for ji = 1:nphy
        plot(tspan,f.Udmsp_phy(ji,:,1),'k')
        plot(tspan,f.Udmsp_phy(ji,:,4),'Color',color(7,:))
        plot(tspan,f.Udmsp_phy(ji,:,3),'Color',color(4,:))
        plot(tspan,f.Udmsp_phy(ji,:,2),'Color',color(2,:))
    end
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP uptake by picophyto [umolS m^{-3} d^{-1}]','fontsize',8)
    grid on
    
%--------------------------------------------------------------------------
% Sulfur fluxes for the article
%--------------------------------------------------------------------------    

figure(fignum(3))
    
    subplot(4,4,1)
    hold on
    plot(tspan,sum(f.leak_DMSP(:,:,1),1),'k-','LineWidth',2)
    plot(tspan,sum(f.mort_to_DMSP(:,:,1),1),'k--','LineWidth',2,'MarkerSize',2)
    plot(tspan,sum(f.graz_to_DMSP(:,:,1),1),'k:','LineWidth',2,'MarkerSize',2)
    legend({'Exudation','Mortality','Grazing  '},'Location','NorthEast')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 10])
    grid on
    
    subplot(4,4,2)
    hold on
    plot(tspan,sum(f.leak_DMSP(:,:,2),1),'-','Color',[0 0.8 0.6],'LineWidth',2)
    plot(tspan,sum(f.mort_to_DMSP(:,:,2),1),'--','Color',[0 0.8 0.6],'LineWidth',2)
    plot(tspan,sum(f.graz_to_DMSP(:,:,2),1),':','Color',[0 0.8 0.6],'LineWidth',2)
    legend({'Exudation','Mortality','Grazing  '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 300])
    set(gca,'ytick',0:50:300)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,3)
    hold on
    plot(tspan,sum(f.leak_DMSP(:,:,3),1),'-','Color',color(4,:),'LineWidth',2)
    plot(tspan,sum(f.mort_to_DMSP(:,:,3),1),'--','Color',color(4,:),'LineWidth',2)
    plot(tspan,sum(f.graz_to_DMSP(:,:,3),1),':','Color',color(4,:),'LineWidth',2)
    legend({'Exudation','Mortality','Grazing  '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 300])
    set(gca,'ytick',0:50:300)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,4)
    hold on
    plot(tspan,sum(f.leak_DMSP(:,:,4),1),'-','Color',color(7,:),'LineWidth',2)
    plot(tspan,sum(f.mort_to_DMSP(:,:,4),1),'--','Color',color(7,:),'LineWidth',2)
    plot(tspan,sum(f.graz_to_DMSP(:,:,4),1),':','Color',color(7,:),'LineWidth',2)
    %legend({['Leakage  ';'Mortality';'Grazing  ']},'Fontsize',11,'Location','NorthWest')
    legend({'Exudation','Mortality','Grazing  '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 300])
    set(gca,'ytick',0:50:300)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,5)
    hold on
    plot(tspan,sum(f.Udmsp_phy(:,:,1),1),'k-','LineWidth',2)
    plot(tspan,f.Udmsp(:,1),'k--','LineWidth',2)
    %legend({['Phytoplankton uptake';'Bacterial uptake    ']},'Fontsize',11,'Location','NorthEast')
    legend({'Phytoplankton uptake','Bacterial uptake    '},'Location','NorthEast')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sinks [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 10])
    grid on
    
    subplot(4,4,6)
    hold on
    plot(tspan,sum(f.Udmsp_phy(:,:,2),1),'-','Color',[0 0.8 0.6],'LineWidth',2)
    plot(tspan,f.Udmsp(:,2),'--','Color',[0 0.8 0.6],'LineWidth',2)
    %legend({['Phytoplankton uptake';'Bacterial uptake    ']},'Fontsize',11,'Location','NorthWest')
    legend({'Phytoplankton uptake','Bacterial uptake    '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sinks [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 300])
    set(gca,'ytick',0:50:300)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,7)
    hold on
    plot(tspan,sum(f.Udmsp_phy(:,:,3),1),'-','Color',color(4,:),'LineWidth',2)
    plot(tspan,f.Udmsp(:,3),'--','Color',color(4,:),'LineWidth',2)
    %legend({['Phytoplankton uptake';'Bacterial uptake    ']},'Fontsize',11,'Location','NorthWest')
    legend({'Phytoplankton uptake','Bacterial uptake    '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sinks [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 300])
    set(gca,'ytick',0:50:300)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,8)
    hold on
    plot(tspan,sum(f.Udmsp_phy(:,:,4),1),'-','Color',color(7,:),'LineWidth',2)
    plot(tspan,f.Udmsp(:,4),'--','Color',color(7,:),'LineWidth',2)
    %legend({['Phytoplankton uptake';'Bacterial uptake    ']},'Fontsize',11,'Location','NorthWest')
    legend({'Phytoplankton uptake','Bacterial uptake    '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP sinks [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 300])
    set(gca,'ytick',0:50:300)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on

    subplot(4,4,9)
    hold on
    plot(tspan,sum(f.leak_DMS(:,:,1),1),'k-','LineWidth',2)
    plot(tspan,f.bac_to_DMS(:,1),'k--','LineWidth',2)
    %legend({['Leakage  ';'Bacteria ']},'Fontsize',11,'Location','NorthEast')
    legend({'Leakage  ','Bacteria '},'Location','NorthEast')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 3])
    grid on
    
    subplot(4,4,10)
    hold on
    plot(tspan,sum(f.leak_DMS(:,:,2),1),'-','Color',[0 0.8 0.6],'LineWidth',2)
    plot(tspan,f.bac_to_DMS(:,2),'--','Color',[0 0.8 0.6],'LineWidth',2)
    %legend({['Leakage  ';'Bacteria ']},'Fontsize',11,'Location','NorthWest')
    legend({'Leakage  ','Bacteria '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 150])
    set(gca,'ytick',0:25:150)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,11)
    hold on
    plot(tspan,sum(f.leak_DMS(:,:,3),1),'-','Color',color(4,:),'LineWidth',2)
    plot(tspan,f.bac_to_DMS(:,3),'--','Color',color(4,:),'LineWidth',2)
    %legend({['Leakage  ';'Bacteria ']},'Fontsize',11,'Location','NorthWest')
    legend({'Leakage  ','Bacteria '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 150])
    set(gca,'ytick',0:25:150)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,12)
    hold on
    plot(tspan,sum(f.leak_DMS(:,:,4),1),'-','Color',color(7,:),'LineWidth',2)
    plot(tspan,f.bac_to_DMS(:,4),'--','Color',color(7,:),'LineWidth',2)
    %legend({['Leakage  ';'Bacteria ']},'Fontsize',11,'Location','NorthWest')
    legend({'Leakage  ','Bacteria '},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sources [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 150])
    set(gca,'ytick',0:25:150)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,13)
    hold on
    plot(tspan,f.DMS_sink(:,1),'k-','LineWidth',2)
    plot(tspan,f.Udms(:,1),'k--','LineWidth',2)
    %legend({['Abiotic sinks   ';'Bacterial uptake']},'Fontsize',11,'Location','NorthEast')
    legend({'Abiotic sinks   ','Bacterial uptake'},'Location','NorthEast')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sinks [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 3])
    grid on
    
    subplot(4,4,14)
    hold on
    plot(tspan,f.DMS_sink(:,2),'-','Color',[0 0.8 0.6],'LineWidth',2)
    plot(tspan,f.Udms(:,2),'--','Color',[0 0.8 0.6],'LineWidth',2)
    %legend({['Abiotic sinks   ';'Bacterial uptake']},'Fontsize',11,'Location','NorthWest')
    legend({'Abiotic sinks   ','Bacterial uptake'},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sinks [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([0 150])
    set(gca,'ytick',0:25:150)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,15)
    hold on
    plot(tspan,f.DMS_sink(:,3),'-','Color',color(4,:),'LineWidth',2)
    plot(tspan,f.Udms(:,3),'--','Color',color(4,:),'LineWidth',2)
    %legend({['Abiotic sinks   ';'Bacterial uptake']},'Fontsize',11,'Location','NorthWest')
    legend({'Abiotic sinks   ','Bacterial uptake'},'Location','NorthWest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sinks [umolS m^{-3} d^{-1}]','fontsize',8) 
    ylim([0 150])
    set(gca,'ytick',0:25:150)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(4,4,16)
    hold on
    plot(tspan,f.DMS_sink(:,4),'-','Color',color(7,:),'LineWidth',2)
    plot(tspan,f.Udms(:,4),'--','Color',color(7,:),'LineWidth',2)
    %legend({['Abiotic sinks   ';'Bacterial uptake']},'Fontsize',11,'Location','NorthWest')
    legend({'Abiotic sinks   ','Bacterial uptake'},'location','northwest')
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS sinks [umolS m^{-3} d^{-1}]','fontsize',8)  
    ylim([0 150])
    set(gca,'ytick',0:25:150)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',4,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
figure(fignum(4))
    
    % Figure of rates: DMSP and DMS turnover rates, DMS yield (DMS prod /
    % DMSP cons)
    
    DMSPd_cons = squeeze(sum(f.Udmsp_phy,1)) + f.Udmsp;
    DMSP_turn = DMSPd_cons ./ DMSPd_mod;
    
    DMS_cons = f.DMS_sink + f.Udms;
    DMS_turn = DMS_cons ./ DMS_mod;
    
    DMS_prod = squeeze(sum(f.leak_DMS + f.graz_to_DMS + f.mort_to_DMS,1)) + f.bac_to_DMS;
    DMSPt_cons = DMSPd_cons + squeeze(sum(f.graz_to_DMSP,1)) * 0.3/0.7; % I add DMSPp assimilation by ZP
    DMS_yield = 100 * DMS_prod ./ DMSPt_cons;
    
    % Net biological production of DMS -> will be negative at the end of
    % the experiment -> "bioloy is a sink of DMS, therefore adding more
    % nutrients willnot increase DMS gas exchange".
    DMS_bio_net = squeeze(sum(f.leak_DMS + f.graz_to_DMS + f.mort_to_DMS,1)) + f.bac_to_DMS - f.Udms;
    
    subplot(2,2,1)
    plot(tspan,DMSP_turn(:,1),'k-','LineWidth',2);
    hold on
    plot(tspan,DMSP_turn(:,4),'-','Color',color(7,:),'LineWidth',2);
    plot(tspan,DMSP_turn(:,3),'-','Color',color(4,:),'LineWidth',2);
    plot(tspan,DMSP_turn(:,2),'-','Color',[0 0.8 0.6],'LineWidth',2);
    xlabel('Time [days]','fontsize',8)
    ylabel('DMSP turnover rate [d^{-1}]','fontsize',8) 
    ylim([0 20])
    set(gca,'ytick',0:2:20)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    %h = legend([h1 h2 h3 h4],{['Control    ';'No-lysis   ';'Late-lysis ';'Early-lysis']},'Fontsize',15,'Location','NorthWest');
    h = legend({'Control    ','No-lysis   ','Late-lysis ','Early-lysis'});
    set(h, 'Position', [0.146 0.8 0.075 0.1])
    grid on
    
    subplot(2,2,2)
    plot(tspan,DMS_turn(:,1),'k-','LineWidth',2)
    hold on
    plot(tspan,DMS_turn(:,4),'-','Color',color(7,:),'LineWidth',2)
    plot(tspan,DMS_turn(:,3),'-','Color',color(4,:),'LineWidth',2)
    plot(tspan,DMS_turn(:,2),'-','Color',[0 0.8 0.6],'LineWidth',2)
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS turnover rate [d^{-1}]','fontsize',8)
    ylim([0 10])
    set(gca,'ytick',0:1:10)
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(2,2,3)
    plot(tspan,DMS_yield(:,1),'k-','LineWidth',2)
    hold on
    plot(tspan,DMS_yield(:,4),'-','Color',color(7,:),'LineWidth',2)
    plot(tspan,DMS_yield(:,3),'-','Color',color(4,:),'LineWidth',2)
    plot(tspan,DMS_yield(:,2),'-','Color',[0 0.8 0.6],'LineWidth',2)
    xlabel('Time [days]','fontsize',8)
    ylabel('DMS yield [%]','fontsize',8) 
    ylim([0 60])
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on
    
    subplot(2,2,4)
    plot(tspan,DMS_bio_net(:,1),'k-','LineWidth',2)
    hold on
    plot(tspan,DMS_bio_net(:,4),'-','Color',color(7,:),'LineWidth',2)
    plot(tspan,DMS_bio_net(:,3),'-','Color',color(4,:),'LineWidth',2)
    plot(tspan,DMS_bio_net(:,2),'-','Color',[0 0.8 0.6],'LineWidth',2)
    xlabel('Time [days]','fontsize',8)
    ylabel('Net biological DMS balance [umolS m^{-3} d^{-1}]','fontsize',8)
    ylim([-80 80])
    yl = ylim;
    ymax = yl(2);
    ymin = yl(1);
    plot([ 8  8],[ymin ymax],'k--','LineWidth',2)
    plot([13 13],[ymin ymax],'k--','LineWidth',2)
    plot([18 18],[ymin ymax],'k--','LineWidth',2)
    plot(0.5:1:7.5,ymin+0.975*(ymax-ymin)*ones(1,8),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    plot(13.5:1:17.5,ymin+0.975*(ymax-ymin)*ones(1,5),'v','MarkerSize',5,'Color',[0.5,0.5,0.5],'MarkerFaceColor',[0.5,0.5,0.5])
    grid on

return

