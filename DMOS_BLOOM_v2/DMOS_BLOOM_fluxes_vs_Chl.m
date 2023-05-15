function [] = DMOS_BLOOM_fluxes_vs_Chl(Chl_obs,DMSPd_obs,DMS_obs,DMSPd_prod,DMS_prod,Ehux_mod,Ehux_death,DMSPd_prod_Ehux,DMS_prod_Ehux,...
                                       color,color4,lineType8,lineType4,fignum)
%DMOS_BLOOM_fluxes_vs_Chl
%   Plot dissolved DMSP and DMS concentration and production as a function
%   of observed chlorophyll concentration

sep = 13; % Separate the first bloom (days 0-13) and the second bloom (days 13-24)

Chl_B1   = nanmean(Chl_obs(:,1:25),2);
Chl_B2   = nanmean(Chl_obs(:,26:end),2);
DMSPd_B1 = nanmean(DMSPd_obs(:,1:25),2);
DMSPd_B2 = nanmean(DMSPd_obs(:,26:end),2);
DMS_B1   = nanmean(DMS_obs(:,1:25),2);
DMS_B2   = nanmean(DMS_obs(:,26:end),2);

Chl_B1_virtual = [Chl_B1(8),0.2*(sum(Chl_B1(1:3))+sum(Chl_B1(5:6))),Chl_B1(4),Chl_B1(7)];
Chl_B2_virtual = [Chl_B2(8),0.2*(sum(Chl_B2(1:3))+sum(Chl_B2(5:6))),Chl_B2(4),Chl_B2(7)];

DMSPd_prod_B1 = mean(DMSPd_prod(1:24*sep,:),1);
DMSPd_prod_B2 = mean(DMSPd_prod(24*sep+1:end,:),1);

DMS_prod_B1 = mean(DMS_prod(1:24*sep,:),1);
DMS_prod_B2 = mean(DMS_prod(24*sep+1:end,:),1);

%Ehux_mod_B1 = mean(Ehux_mod(1:24*sep,:),1);
Ehux_mod_B2 = mean(Ehux_mod(24*sep+1:end,:),1);

%Ehux_death_B1 = mean(Ehux_death(1:24*sep,:),1);
Ehux_death_B2 = mean(Ehux_death(24*sep+1:end,:),1);

%DMSPd_prod_Ehux_B1 = mean(DMSPd_prod_Ehux(1:24*sep,:),1);
DMSPd_prod_Ehux_B2 = mean(DMSPd_prod_Ehux(24*sep+1:end,:),1);

%DMS_prod_Ehux_B1 = mean(DMS_prod_Ehux(1:24*sep,:),1);
DMS_prod_Ehux_B2 = mean(DMS_prod_Ehux(24*sep+1:end,:),1);

ref_slope_DMSPd = mean(DMSPd_B1([1,2,3,5,6])) / mean(Chl_B1([1,2,3,5,6]));
ref_slope_DMS   = mean(DMS_B1([1,2,3,5,6]))   / mean(Chl_B1([1,2,3,5,6]));
ref_slope_DMSPd_prod = DMSPd_prod_B1(2) / Chl_B1_virtual(2);
ref_slope_DMS_prod   = DMS_prod_B1(2)   / Chl_B1_virtual(2);

figure(fignum(1))
subplot(2,2,1)
hold on
for ji=1:8
    plot(Chl_B1(ji),DMSPd_B1(ji),lineType8{ji},'Markersize',12,'Color',color(ji,:),'LineWidth',2)
    plot(Chl_B2(ji),DMSPd_B2(ji),lineType8{ji},'Markersize',12,'Color',color(ji,:),'MarkerFaceColor',color(ji,:))
end
plot([0 16],ref_slope_DMSPd*[0 16],'k--','LineWidth',2)
xlabel('Average observed chlorophyll concentration [mg m^{-3}]')
ylabel('Average observed DMSPd concentration [umol m^{-3}]')
xlim([0 16])
ylim([0 50])
grid on
subplot(2,2,2)
hold on
for ji = 1:8
    plot(Chl_B1(ji),DMS_B1(ji),lineType8{ji},'Markersize',12,'Color',color(ji,:),'LineWidth',2)
    plot(Chl_B2(ji),DMS_B2(ji),lineType8{ji},'Markersize',12,'Color',color(ji,:),'MarkerFaceColor',color(ji,:))
end
plot([0 16],ref_slope_DMS*[0 16],'k--','LineWidth',2)
xlabel('Average observed chlorophyll concentration [mg m^{-3}]')
ylabel('Average observed DMS concentration [umol m^{-3}]')
xlim([0 16])
ylim([0 50])
grid on
subplot(2,2,3)
hold on
for ji = 1:4
    plot(Chl_B1_virtual(ji),DMSPd_prod_B1(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'LineWidth',2)
    plot(Chl_B2_virtual(ji),DMSPd_prod_B2(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'MarkerFaceColor',color4(ji,:))
end
plot([0 16],ref_slope_DMSPd_prod*[0 16],'k--','LineWidth',2)
xlabel('Average observed chlorophyll concentration [mg m^{-3}]')
ylabel('Average model DMSPd production [nM d^{-1}]')
xlim([0 16])
ylim([0 120])
grid on
subplot(2,2,4)
hold on
for ji = 1:4
    plot(Chl_B1_virtual(ji),DMS_prod_B1(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'LineWidth',2)
    plot(Chl_B2_virtual(ji),DMS_prod_B2(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'MarkerFaceColor',color4(ji,:))
end
plot([0 16],ref_slope_DMS_prod*[0 16],'k--','LineWidth',2)
xlabel('Average observed chlorophyll concentration [mg m^{-3}]')
ylabel('Average model DMS production [nM d^{-1}]')
xlim([0 16])
ylim([0 50])
grid on

figure(fignum(2))
subplot(2,2,1)
hold on
for ji = 1:4
    plot(Ehux_mod_B2(ji),DMSPd_prod_Ehux_B2(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'MarkerFaceColor',color4(ji,:))
end
xlabel('Average model E. huxleyi biomass [uM]')
ylabel('Average model Ehux DMSPd production [nM d^{-1}]')
grid on
subplot(2,2,2)
hold on
for ji = 1:4
    plot(Ehux_mod_B2(ji),DMS_prod_Ehux_B2(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'MarkerFaceColor',color4(ji,:))
end
xlabel('Average model E. huxleyi biomass [uM]')
ylabel('Average model Ehux DMS production [nM d^{-1}]')
grid on
subplot(2,2,3)
hold on
for ji = 1:4
    plot(Ehux_death_B2(ji),DMSPd_prod_Ehux_B2(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'MarkerFaceColor',color4(ji,:))
end
xlabel('Average model E. huxleyi death flux [uM d^{-1}]')
ylabel('Average model Ehux DMSPd production [nM d^{-1}]')
grid on
subplot(2,2,4)
hold on
for ji = 1:4
    plot(Ehux_death_B2(ji),DMS_prod_Ehux_B2(ji),lineType4{ji},'Markersize',12,'Color',color4(ji,:),'MarkerFaceColor',color4(ji,:))
end
xlabel('Average model E. huxleyi death flux [uM d^{-1}]')
ylabel('Average model Ehux DMS production [nM d^{-1}]')
grid on

%dbstop

end

