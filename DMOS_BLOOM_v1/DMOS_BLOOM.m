%********************************************************************
% SUMMIT (Novel roles of dimethylated sulphur in marine microbial
% food webs) version 1.0
% Code written by Guillaume Le Gland and Sergio M. Vallina
%********************************************************************

%********************************************************************
% PROGRAM: SUMMIT-BERGEN.M
%********************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FILE HEADER -- LOAD PACKAGES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************************************************
more off
close all
clear all
format short g
%--------------------------------------------------------------------
addpath(genpath('TOOLBOX'));
%--------------------------------------------------------------------
%====================================================================
%--------------------------------------------------------------------
[mypackages] = myheadloadpackages; %Structure array to pass on my head pkg as input argument to functions.
%--------------------------------------------------------------------
%MY PACKAGES FOR PLOTING:
subplot_funhan  = mypackages.subplot;
colorbar_funhan = mypackages.colorbar;
verticales = mypackages.verticales;
horizontal = mypackages.horizontal;
%--------------------------------------------------------------------
%JUST CHECKING IF PLOTING WORKS OKAY:
A256 = peaks(256);
A128x256 = A256(1:2:256,:);
A = A128x256;
Amin = min(A(:));
Amax = max(A(:));
fignum = 1001;
[hcbar] = DMOS_BLOOM_subplotesting(A,Amin,Amax,fignum,mypackages);
pause(0.5)
close all
%--------------------------------------------------------------------
%====================================================================
%********************************************************************
global t0 deltat ndays nyear tmax tspan 
%...................................................................................
global Jpico Jnano Jehux Jzoo Jbac Jdin Jon
global Jdmspd Jdms Jfdc
global keyGrowthDMSP keyGrowthDMS keyVarDMSUptake
global day Ehux_obs nano_obs pico_obs Syn_obs bac_obs EhV_obs MpV_obs cil_obs thr_obs NO3_obs Chl_obs DMS_obs DMSPd_obs DMSPp_obs
%...................................................................................
global FLUXES % Sources and sinks outputs

tic 

[keyInvPlankton,keyInvSulfur,keyObsTolerance,keyGrowthDMSP,keyGrowthDMS,keyVarDMSUptake] = DMOS_BLOOM_keys(); 

[param0,paramOptX,paramOptY,X0,Y0,xtl,ytl,Xtol,Ytol,lowBndX,uppBndX,lowBndY,uppBndY] = DMOS_BLOOM_parameters();

nx = length(X0);
ny = length(Y0);

%%%%%%%%%%%%%%%%%%%%%
%TEMPORAL RESOLUTION:
%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................
nyear = 1;
deltat = param0.deltat;
ndays = param0.ndays;
t0 = deltat;
tmax = ndays*nyear; 
tspan = t0:deltat:tmax; 
%........................................................................
nsteps = length(tspan); % Should be equal to (ndays/deltat)*nyear
ttmax = ndays; %[days]
%........................................................................
%========================================================================
%%%%%%%%%%%%%%%%%%%
%DATA READING
%%%%%%%%%%%%%%%%%%%

[day,Ehux_obs,nano_obs,pico_obs,Syn_obs,bac_obs,EhV_obs,MpV_obs,cil_obs,thr_obs,IEC_obs,...
 NO3_obs,Chl_obs,DMS_obs,DMSPd_raw,DMSPd_obs,DMSPp_obs] = DMOS_BLOOM_load_data();

fignum = [10,20,30,40,50,60,70];
DMOS_BLOOM_plot_data(day,Ehux_obs,nano_obs,pico_obs,Syn_obs,bac_obs,EhV_obs,...
    cil_obs,thr_obs,IEC_obs,NO3_obs,Chl_obs,DMS_obs,DMSPd_raw,DMSPd_obs,DMSPp_obs,fignum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINE THE ROWS CORRESPONDING TO EACH STATE-VARIABLE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%FOR DISCRETE TRAIT MODEL:
%------------------------------------------------------------------------
%[NOTE: Jphy, Jzoo, Jdin, Jpom]
%------------------------------------------------------------------------
%........................................................................
nvar = [0,1,2,1,5,1,1,1,1,1,1];
%........................................................................
varname{1} = 'pico';
varname{2} = 'nano';
varname{3} = 'ehux';
varname{4} = 'zoo';
varname{5} = 'bac';
varname{6} = 'din';
varname{7} = 'on';
varname{8} = 'dms';
varname{9} = 'dmspd';
varname{10} = 'fdc'; % Fraction of DMS consumers
%........................................................................
nsum1 = 0;
nsum2 = 0;
%........................................................................
for j = 1:length(nvar)-1
    %....................................................................
    nsum1 = nsum1 + nvar(j);
    nsum2 = nsum2 + nvar(j+1);
    %....................................................................
    %Jvarj = (ndepths*(nsum1))+1:(ndepths*(nsum2));
    Jvarj = (nsum1)+1:(nsum2);
    %....................................................................
    Jstr = ['J',varname{j}];
    %....................................................................
    %% assignin('base',Jstr,Jvarj); %DO NOT USE THIS OLD APPROACH
    %....................................................................
    myassign = [Jstr,' = Jvarj;']; %BETTER TO USE THIS NEW APPROACH FOR GLOBAL VARIABLES
    %....................................................................
    eval(myassign)
    %....................................................................
end

%--------------------------------------------------------------------------
% Optimization ("inverse modeling")
% Takes a few tens of minutes in Matllab but several hours in Octave
%--------------------------------------------------------------------------

%..........................................................................
% 4 simulations:
% 1) Control without DIN input -> Fjord
% 2) Standard with DIN input and bubbling -> Bags 1, 2, 3, 5 and 6
% 3) DIN input + bubbling + extra mortality on E.hux from day 17 on -> bag 4
% 4) DIN input + bubbling + extra mortality on E.hux from day 11 on -> bag 7
%..........................................................................

bag_exp = [2,2,2,0,2,2,0,1];

raw_obs_plankton = cat(3,pico_obs,nano_obs,Ehux_obs,cil_obs,bac_obs,NO3_obs);

% Account for NH4 and NO2
offset = 2*nanmean(NO3_obs(8,:)); 
raw_obs_plankton(:,:,6)   = raw_obs_plankton(:,:,6) + offset;
weight_obs_plankton = [1,1,1,0,2,1];

raw_obs_sulfur = cat(3,DMS_obs,DMSPd_obs,DMSPp_obs);
weight_obs_sulfur = [1,1,1]; % Weight of each variable in the sulfur cost function

[obsPlankton,obsPlanktonTol] = DMOS_BLOOM_format_data(raw_obs_plankton,weight_obs_plankton,bag_exp,keyObsTolerance);
[obsSulfur  ,obsSulfurTol  ] = DMOS_BLOOM_format_data(raw_obs_sulfur  ,weight_obs_sulfur  ,bag_exp,keyObsTolerance);

day_obs = day(1,:);

options = optimset('MaxIter',2000);

func_plankton = @(x) DMOS_BLOOM_cost_plankton(tspan,x,X0,Xtol,...
      param0,paramOptX,day_obs,obsPlankton,obsPlanktonTol);

%..........................................................................
% Optimization of the plankton model
%..........................................................................

if strcmp(keyInvPlankton,'yes')
    [X,RESNORMX,RESIDUALX,EXITFLAGX,OUTPUTX,LAMBDAX,JACOBIANX] = lsqnonlin(func_plankton,X0,lowBndX,uppBndX,options);
    % Estimate of uncertainty (following Le Gland et al., 2017)
    % Mean square of parameters
    etax = mean(RESIDUALX(:).^2);
    % Factor involving the Jacobian matrix
    sensx = inv(JACOBIANX' * JACOBIANX);
    % Numerical factor acccounting for the number of data and parameters
    ndp = length(RESIDUALX); 
    facx = ndp / (ndp - nx);
    varxx = facx * sensx * etax;
    varxx_norm = zeros(size(varxx));
    for ji = 1:nx
        for jj = 1:nx
            varxx_norm(ji,jj) = varxx(ji,jj) / sqrt(varxx(ji,ji)*varxx(jj,jj));
        end
    end
    varxx_diag = diag(varxx);
    stdx = sqrt(varxx_diag);
    % Plot parameters with their error bars
    figure(100)
    %opt = {'g','Linew',2,'LineS','none'};
    %errorbar(1:nx,X,stdx,stdx,opt{:});
    errorbar(1:nx,X,stdx,stdx,'g.')
    hold on
    plot(1:nx,zeros(1,nx),'k')
    xlabel('Plankton parameters','fontsize',15)
    set(gca,'xtick',1:nx);
    set(gca,'xticklabel',xtl,'fontsize',10);
    ylabel('Confidence interval','fontsize',15)
    % Plot correlations
    figure(110)
    imagesc(varxx_norm,[-1,+1])
    colorbar
    xlabel('Plankton parameters','fontsize',15)
    set(gca,'xtick',1:nx);
    set(gca,'xticklabel',xtl,'fontsize',10);
    ylabel('Plankton parameters','fontsize',15)
    set(gca,'ytick',1:nx);
    set(gca,'yticklabel',xtl,'fontsize',10);
    title('Plankton correlation matrix','fontsize',20)
    colorbar
else
    % Results of standard optimization. Can be used directly.
    X = [0.662 0.673 0.602 0.654 3.52 8.32 0.587 4.32 2.34 2.69 1.81 1.81 2.21 2.2 2.73 2.46 3.01 1.15 0.938 0.598 4.27e-010 0.557];
end
% Update param0 for later sulfur optimization
fnVarX = fieldnames(paramOptX);
nVarX = numel(fnVarX);
disp(['X:  ', mat2str(X',3)])
for ji=1:nVarX
    arr = paramOptX.(fnVarX{ji});
    sz = length(arr);
    for n=1:sz
        if arr(n) > 0
            param0.(fnVarX{ji})(n) = X(arr(n));
        end
    end
end
%param0.Kon = param0.mub0 /  param0.Kon;
%param0.Kz  = param0.muz0 ./ param0.Kz;

%..........................................................................
% Optimization of the DMSP & DMS model
% It is not used in the article and I do not recommend its use, due to the
% low amount of data and large number of parameters.
% For instance, the algorithm tries to remove DMS uptake by co-metabolizers
% ("generalists"), and imposes high rates of DMSP leakage.
%..........................................................................

func_sulfur = @(y) DMOS_BLOOM_cost_sulfur(tspan,y,Y0,Ytol,...
      param0,paramOptY,day_obs,obsSulfur,obsSulfurTol);

if strcmp(keyInvSulfur,'yes')
    [Y,RESNORMY,RESIDUALY,EXITFLAGY,OUTPUTY,LAMBDAY,JACOBIANY] = lsqnonlin(func_sulfur,Y0,lowBndY,uppBndY);
    etay = mean(RESIDUALY(:).^2);
    sensy = inv(JACOBIANY' * JACOBIANY);
    %sensy = JACOBIANY' * JACOBIANY;
    nds = length(RESIDUALY); 
    facy = nds / (nds - ny);
    varyy = facy * sensy * etay;
    varyy_norm = zeros(size(varyy));
    for ji = 1:ny
        for jj = 1:ny
            varyy_norm(ji,jj) = varyy(ji,jj) / sqrt(varyy(ji,ji)*varyy(jj,jj));
        end
    end
    varyy_diag = diag(varyy);
    stdy = sqrt(varyy_diag);
    figure(120)
    %opt = {'r','Linew',2,'LineS','none'};
    %errorbar(1:ny,100*Y./Y0,100*stdy./Y0,100*stdy./Y0,opt{:})
    errorbar(1:ny,Y,stdy,stdy,'r.')
    xlabel('Sulfur parameters','fontsize',15)
    set(gca,'xtick',1:ny);
    set(gca,'xticklabel',ytl,'fontsize',10);
    ylabel('Confidence interval (% of first guess)','fontsize',15)
    hold on
    plot(1:ny,zeros(1,ny),'k')
    figure(130)
    imagesc(varyy_norm,[-1,+1])
    colorbar
    xlabel('Sulfur parameters','fontsize',15)
    set(gca,'xtick',1:ny);
    set(gca,'xticklabel',ytl,'fontsize',10);
    ylabel('Sulfur parameters','fontsize',15)
    set(gca,'ytick',1:ny);
    set(gca,'yticklabel',ytl,'fontsize',10);
    title('Sulfur correlation matrix','fontsize',20)
    colorbar
    % Update param0 for later sulfur optimization
    fnVarY = fieldnames(paramOptY);
    nVarY = numel(fnVarY);
    for ji=1:nVarY
        arr = paramOptY.(fnVarY{ji});
        sz = length(arr);
        for n=1:sz
            if arr(n) > 0
                param0.(fnVarY{ji})(n) = Y(arr(n));
            end
        end
    end
    param0.mort_DMSP_frac = 1 - param0.mort_DMS_frac;
    graz_frac = 0.7;
    param0.graz_DMS_frac = graz_frac * param0.mort_DMS_frac;
    param0.gras_DMSP_frac = graz_frac * (1 - param0.graz_DMS_frac);
    disp(['Y:  ', mat2str(Y',3)])
else
        Y = Y0;
end

%...............................................................................

param0.mort_DMSP_frac = 1 - param0.mort_DMS_frac;
graz_frac = 0.7;
param0.graz_DMS_frac = graz_frac * param0.mort_DMS_frac;
param0.graz_DMSP_frac = graz_frac - param0.graz_DMS_frac;

Vdisc0 = [param0.pico0;
          param0.nano0;
          param0.Ehux0;
          param0.zoo0;
          param0.bac0;
          param0.DIN0;
          param0.ON0;
          param0.dms0;
          param0.dmspd0;
          param0.fdc0];

% Remove fields that the model does not need
param0 = rmfield(param0,{'pico0','nano0','Ehux0','zoo0','bac0',...
    'DIN0','dmspd0','dms0'});

nout = length(Vdisc0);
Vout = zeros(nsteps,nout,4);

% Output the fluxes (global variables)
param0.outputFluxes = 'yes';

%--------------------------------------------------------------------------
% Model outputs
%--------------------------------------------------------------------------

if strcmp(param0.outputFluxes,'yes')
    nphy = 4;
    nzoo = 5;
    nsim = 4; % Number of virtual bags
    FLUXES.uptake_DIN_phy = zeros(nphy,nsteps,nsim);
    FLUXES.mort_phy       = zeros(nphy,nsteps,nsim);
    FLUXES.ON_to_bac      = zeros(nsteps,nsim);
    FLUXES.remin_ON       = zeros(nsteps,nsim);
    FLUXES.mort_bac       = zeros(nsteps,nsim);
    FLUXES.graz_zoo       = zeros(nzoo,nsteps,nsim);
    FLUXES.prey_to_zoo    = zeros(nzoo,nsteps,nsim);
    FLUXES.mort_zoo       = zeros(nzoo,nsteps,nsim);
    FLUXES.leak           = zeros(nzoo,nsteps,nsim);
    FLUXES.sink_ON        = zeros(nsteps,nsim);
    FLUXES.leak_DMS       = zeros(nphy,nsteps,nsim);
    FLUXES.bac_to_DMS     = zeros(nsteps,nsim);
    FLUXES.mort_to_DMS    = zeros(nphy,nsteps,nsim);
    FLUXES.graz_to_DMS    = zeros(nphy,nsim);
    FLUXES.leak_DMSP      = zeros(nphy,nsteps,nsim);
    FLUXES.mort_to_DMSP   = zeros(nphy,nsteps,nsim);
    FLUXES.graz_to_DMSP   = zeros(nphy,nsteps,nsim);
    FLUXES.DMS_sink       = zeros(nsteps,nsim);
    FLUXES.Udms           = zeros(nsteps,nsim);
    FLUXES.Udmsp          = zeros(nsteps,nsim);
    FLUXES.Udmsp_phy      = zeros(nphy,nsteps,nsim);
    % FLUXES.DMSPp          = zeros(nsteps,nsim);
end

% 1) Control without DIN input -> Fjord
param0.inputCode = 0;
param0.lysisTime = -1;
param0.isim = 1;
Vout(:,:,1) = ode4(@DMOS_BLOOM_ode45eqs,tspan,Vdisc0,param0);

% 2) Standard with DIN input days 0-8 and 13-18 -> Bags 1, 2, 3, 5 and 6
param0.inputCode = 1;
param0.lysisTime = -1;
param0.isim = 2;
Vout(:,:,2) = ode4(@DMOS_BLOOM_ode45eqs,tspan,Vdisc0,param0); % 21

% 3) DIN input + extra mortality on E.hux from day 17 on (equivalent of bag 4)
param0.inputCode = 1;
param0.lysisTime = 17;
param0.isim = 3;
Vout(:,:,3) = ode4(@DMOS_BLOOM_ode45eqs,tspan,Vdisc0,param0);

% 4) DIN input + extra mortality on E.hux from day 11 on (equivalent of bag 7)
param0.inputCode = 1;
param0.lysisTime = 11;
param0.isim = 4;
Vout(:,:,4) = ode4(@DMOS_BLOOM_ode45eqs,tspan,Vdisc0,param0);

% Plankton outputs
pico_mod = squeeze(Vout(:,Jpico,:));
nano_mod = squeeze(sum(Vout(:,Jnano,:),2));
nano1_mod = squeeze(Vout(:,Jnano(1),:));
nano2_mod = squeeze(Vout(:,Jnano(2),:));
Ehux_mod = squeeze(Vout(:,Jehux,:));
zoo_mod  = squeeze(sum(Vout(:,Jzoo,:),2));
bac_mod  = squeeze(Vout(:,Jbac,:));
DIN_mod  = squeeze(Vout(:,Jdin,:));
ON_mod   = squeeze(sum(Vout(:,Jon,:),2));

% Sulfur outputs
DMS_mod   = squeeze(Vout(:,Jdms,:));
DMSPd_mod = squeeze(Vout(:,Jdmspd,:));
FDC_mod    = squeeze(Vout(:,Jfdc,:));
DMSPp_mod = zeros(size(DIN_mod));
for ji = 1:size(DIN_mod,1)
    for jj = 1:size(DIN_mod,2)
        lim_DIN = DIN_mod(ji,jj) ./ (DIN_mod(ji,jj) + param0.Kn); % vector
        DMSP_ratio = param0.thetaSN_P .* (1-(1-param0.stress_fac).*lim_DIN); % vector
        DMSP_pico  = DMSP_ratio(1) * (1000*pico_mod(ji,jj));
        DMSP_nano1 = DMSP_ratio(2) * (1000*nano1_mod(ji,jj));
        DMSP_nano2 = DMSP_ratio(3) * (1000*nano2_mod(ji,jj));
        DMSP_nano  = DMSP_nano1 + DMSP_nano2;
        DMSP_Ehux  = DMSP_ratio(4) * (1000*Ehux_mod(ji,jj));
        % DMSPp_mod(ji,jj)  = FLUXES.DMSPp(ji,jj);
        DMSPp_mod(ji,jj) = DMSP_pico + DMSP_nano + DMSP_Ehux;
    end
end

ChlCrat = 0.026; % mgChl / mgC
ChlNrat = ChlCrat * 6.625;
Chl_mod = (pico_mod + nano_mod + Ehux_mod) * 12.0 * ChlNrat;  

% Correlation coefficient and Lin's CCC for each variable (mod vs obs, No-lysis virtual bag)
mask = ~isnan(obsPlankton(2,:,5));
mod = bac_mod(max(1,24*day(1,mask)),2);
obs = obsPlankton(2,mask,5)';
corr_bac = corr(mod,obs);
CCC_bac  = corr_bac * (2 * std(mod) * std(obs)) / (var(mod) + var(obs) + (mean(mod)-mean(obs))^2);
%..........................................................................
mask = ~isnan(obsPlankton(2,:,1));
mod = pico_mod(max(1,24*day(1,mask)),2);
obs = obsPlankton(2,mask,1)';
corr_pico = corr(mod,obs);
CCC_pico  = corr_pico * (2 * std(mod) * std(obs)) / (var(mod) + var(obs) + (mean(mod)-mean(obs))^2);
%..........................................................................
mask = ~isnan(obsPlankton(2,:,2));
mod = nano_mod(max(1,24*day(1,mask)),2);
obs = obsPlankton(2,mask,2)';
corr_nano = corr(mod,obs);
CCC_nano  = corr_nano * (2 * std(mod) * std(obs)) / (var(mod) + var(obs) + (mean(mod)-mean(obs))^2);
%..........................................................................
mask = ~isnan(obsPlankton(2,:,3));
mod = Ehux_mod(max(1,24*day(1,mask)),2);
obs = obsPlankton(2,mask,3)';
corr_Ehux = corr(mod,obs);
CCC_Ehux  = corr_Ehux * (2 * std(mod) * std(obs)) / (var(mod) + var(obs) + (mean(mod)-mean(obs))^2);
%..........................................................................
mask = ~isnan(obsPlankton(2,:,6));
mod = DIN_mod(max(1,24*day(1,mask)),2);
obs = obsPlankton(2,mask,6)';
corr_DIN = corr(mod,obs);
CCC_DIN  = corr_nano * (2 * std(mod) * std(obs)) / (var(mod) + var(obs) + (mean(mod)-mean(obs))^2);
%..........................................................................

color = [0.5 0 0; 0.7 0.7 0; 0.7 0 0.7; 1 0 0; 0 0 0.5; 0 0.7 0.7; 0 0 1; 0 0 0];
color4 = [0 0 0; 0 0.8 0.6; 1 0 0; 0 0 1];
faceColor = [1 1 1; 1 1 1; 1 1 1; 1 0 0; 1 1 1; 1 1 1; 0 0 1; 0 0 0];

lineType  = {'o--';'s--';'d--';'o- ';'v--';'^--';'s- ';'d- '};
lineType8 = {'o';'s';'d';'o';'v';'^';'s';'d'};
lineType4 = {'d';'p';'o';'s'};

DMSPd_prod = squeeze(sum(FLUXES.leak_DMSP + FLUXES.graz_to_DMSP + FLUXES.mort_to_DMSP,1));
DMS_prod = squeeze(sum(FLUXES.leak_DMS + FLUXES.graz_to_DMS + FLUXES.mort_to_DMS,1)) + FLUXES.bac_to_DMS;

Ehux_death = squeeze(FLUXES.mort_phy(4,:,:)+FLUXES.graz_zoo(5,:,:));

DMSPd_prod_Ehux = squeeze(FLUXES.leak_DMSP(4,:,:) + FLUXES.graz_to_DMSP(4,:,:) + FLUXES.mort_to_DMSP(4,:,:));
DMS_prod_Ehux   = squeeze(FLUXES.leak_DMS(4,:,:)  + FLUXES.graz_to_DMS(4,:,:)  + FLUXES.mort_to_DMS(4,:,:)) + param0.alpha2 * DMSPd_prod_Ehux;

fignum = [200,210,220];
DMOS_BLOOM_plot_model(zoo_mod,Chl_mod,pico_mod,nano_mod,nano1_mod,nano2_mod,Ehux_mod,bac_mod,...
                      DIN_mod,ON_mod,DMSPp_mod,DMSPd_mod,DMS_mod,FDC_mod,...
                      obsPlankton,cil_obs,Chl_obs,pico_obs,nano_obs,Ehux_obs,bac_obs,NO3_obs,...
                      obsSulfur,DMSPp_obs,DMSPd_obs,DMS_obs,offset,tspan,day,color,fignum)

if strcmp(param0.outputFluxes,'yes')
    fignum = [300,310,320,330];
    DMOS_BLOOM_plot_fluxes(FLUXES,DMSPd_mod,DMS_mod,tspan,color,fignum)
end

fignum = [400,410];
DMOS_BLOOM_fluxes_vs_Chl(Chl_obs,DMSPd_obs,DMS_obs,DMSPd_prod,DMS_prod,Ehux_mod,Ehux_death,DMSPd_prod_Ehux,DMS_prod_Ehux,...
                         color,color4,lineType8,lineType4,fignum)

[day_Meth,Meth_free,Meth_part] = DMOS_BLOOM_load_DMS_consumers();

fignum = 420;
DMOS_BLOOM_plot_DMS_consumers(day_Meth,Meth_free,Meth_part,ndays,color,faceColor,lineType,fignum)

return