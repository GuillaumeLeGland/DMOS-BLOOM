function [] = DMOS_BLOOM_plot_PAR(PAR,Chl_obs,day_obs,ndays,color,fignum)
% Plot surface and average PAR in the mesocosm bags over the duration of
% the experiment

% Convert chlorophyll observations to a 24-day array
nbags = size(Chl_obs,1);
ndat  = size(Chl_obs,2);
Chl24 = zeros(nbags,ndays);
for bag = 1:nbags
    day = 1;
    nobs = 0; % Number of observations at day d
    for dat = 1:ndat
        if day_obs(dat) >= day && day ~= ndays
            Chl24(bag,day) = Chl24(bag,day) / max(nobs,1);
            day = day + 1;
            nobs = 0;
        end
        if ~isnan(Chl_obs(bag,dat))
            Chl24(bag,day) = Chl24(bag,day) + Chl_obs(bag,dat);
            nobs = nobs + 1;
        end       
    end
end

% Estimate PAR averaged over a 4-m column 
PAR_ave = zeros(nbags,ndays);
kw = 0.04; % PAR attenuation due to water [1/m], same as in Le Gland et al. (2021) and Vallina et al. (2008)
kchl = 0.03; %0.03; % PAR attenuation due to chlorophyll [m²/mgChl] (upper estimate from Kirk 2011 Fig. 3.11 for natural phytoplankton communities)
zmax = 4; % Maximum depth [m]
for bag = 1:nbags
    for day = 1:ndays
        Kd = kw + kchl*Chl24(bag,day);
        PAR_ave(bag,day) = PAR(day) * (1/(zmax*Kd)) * (1 - exp(-zmax*Kd));
    end
end

ymax = 160;

PAR_ave_4 = zeros(nbags,ndays);
PAR_ave_4(4,:) = PAR_ave(4,:);
PAR_ave_7 = zeros(nbags,ndays);
PAR_ave_7(7,:) = PAR_ave(7,:);

figure(fignum)
hold on
bar(1:ndays,PAR',0.8,'facecolor',[0.9 0.9 0])
bar(1:ndays,PAR_ave(1:7,:)',0.7,'facecolor',[0.3 0.3 0.3])
bar(1:ndays,PAR_ave_4(1:7,:)',0.7,'facecolor',[1 0 0])
bar(1:ndays,PAR_ave_7(1:7,:)',0.7,'facecolor',[0 0 1])

ylim([0 ymax])
plot([ 8.5  8.5],[0 ymax],'k--','LineWidth',5)
plot([13.5 13.5],[0 ymax],'k--','LineWidth',5)
plot([18.5 18.5],[0 ymax],'k--','LineWidth',5)
plot(1:1:8,0.975*ymax*ones(1,8),'v','MarkerSize',12,'Color',[0.3,0.3,0.3],'MarkerFaceColor',[0.3,0.3,0.3])
plot(14:1:18,0.975*ymax*ones(1,5),'v','MarkerSize',12,'Color',[0.3,0.3,0.3],'MarkerFaceColor',[0.3,0.3,0.3])
xlabel('Time [days]','fontsize',8)
ylabel('Surface PAR','fontsize',8)
grid on

return