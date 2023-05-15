function [Vdot] = DMOS_BLOOM_ode45eqs(iTime,V0,varargin)
global icounter 
global deltat
global Jpico Jnano Jehux Jzoo Jbac Jdin Jon 
global Jdmspd Jdms Jfdc
global jday
global keyGrowthDMSP keyGrowthDMS keyVarDMSUptake keyTemp keyPAR
%...................................................................................

%%%%%%%%%%%%%%%%%
%STATE VARIABLES:
%%%%%%%%%%%%%%%%%
%===================================================================================
%...................................................................................
pico   = V0(Jpico);
nano   = V0(Jnano);
Ehux   = V0(Jehux);
zoo    = V0(Jzoo);
bac    = V0(Jbac);
DIN    = V0(Jdin);
ON     = V0(Jon);
DMS    = V0(Jdms);
DMSPd  = V0(Jdmspd);
FDC     = V0(Jfdc); % Fraction of DMS-consuming heterotrophic bacteria
%...................................................................................
DIN   = max(0,DIN);
ON    = max(0,ON);
DMSPd = max(0,DMSPd);
DMS   = max(0,DMS);
%...................................................................................

% Load parameters
p = varargin{1};
deltat    = p.deltat;

% Other parameters
input_DIN = zeros(24,1); % DIN input to the bags
input_ON  = zeros(24,1); % ON input in the fjord
if p.inputCode == 0
    input_ON(:,end)  = sum(p.ON0) * p.wON;
elseif p.inputCode > 0
    input_DIN(1:8)   = input_DIN(1:8) + 1.6;
    input_DIN(14:18) = input_DIN(14:18) + 1.6;
end

%%%%%%%%%%%%%%%%%%%
%DAY OF SIMULATION:
%%%%%%%%%%%%%%%%%%%
%[jday,newday] = DMOS_BLOOM_daycounter(jday,iTime);
jday = DMOS_BLOOM_daycounter(jday,iTime);
t = ceil(iTime);

%..........................................................................
% Temperature effect on metabolic rates (Q10 approach)
if strcmp(keyTemp,'yes')
    temp = p.temp(p.inputCode+1,jday);
    fT_a = exp(0.056*(temp-13.0)); % Autotrophic processes
    fT_h = exp(0.092*(temp-13.0)); % heterotrophic processes
else
    % If temperature effect is not taken into account
    fT_a = 1.0;
    fT_h = 1.0;
end

% Light effect on metabolic rates
if strcmp(keyPAR,'yes')
    PAR  = p.PAR(jday);
    ChlCrat = 0.026; % mgChl / mgC
    ChlNrat = ChlCrat * 6.625;
    Chl     = (pico + sum(nano(:)) + Ehux) * 12.0 * ChlNrat;
    Kd = 0.04 + 0.03*Chl; % Downwelling light attenuation
    % Average light over the depth of the mesocosm (4 meters)
    jPAR = PAR * (1/(4*Kd)) * (1 - exp(-4*Kd)); 
    % Light limitation (Follows et al. 2007 + Le Gland et al. 2021)
    InhFac = 12;
    Isat = 60; % From Vallina et al. (2008)
    Kpar   = log(InhFac+1) / Isat;
    Kinhib = Kpar / InhFac;
    Fmax   = (Kpar+Kinhib)/Kpar * exp( -(Kinhib/Kpar) * log(Kinhib/(Kpar+Kinhib)) );  
    % Qpar is normalized to have a maximum of 1 at Isat
    Qpar = Fmax * (1 - exp(-Kpar*jPAR)) .* exp(-Kinhib*jPAR);
    % UV multiplicative factor on DMSP and DMS release
    UV_fac = jPAR/60;
    % Photoinhibition of bacteria
    UV_bac = 1 - (1-p.stress_fac)*(jPAR/120);
else
    % If light effect is not taken into account
    Qpar = 1.0;
    UV_fac = 1.0;
    UV_bac = 1.0;
end
%..........................................................................

%........................................................................
%Ineg = find(VNPZD0 < 0); 
%...................................................................................
%===================================================================================

%%%%%%%%%%%
%SHOW TIME:
%%%%%%%%%%%
%===================================================================================
%...................................................................................
icounter = floor(iTime/p.deltat);
%........................................................................
%===================================================================================

%..........................................................................
% Plankton model
%..........................................................................

phy = [pico;nano;Ehux];

uptake_DIN_phy = p.mup0 .* phy .* DIN ./ (DIN + p.Kn) .* (fT_a.*Qpar);
uptake_DIN_phy_tot = sum(uptake_DIN_phy(:));

lim_DIN = DIN ./ (DIN + p.Kn );
lim_ON  = ON  ./ (ON  + p.Kon);

prey = [bac;phy];
% One zooplankton for each prey type
grazZoo   = (fT_h) .* p.muz0 .* zoo .* prey ./ (p.Kz + prey); % Holling type II
exudZoo   = (1 - p.betaz) .* grazZoo;
preyToZoo =      p.betaz   .*grazZoo;  
exudToDIN = sum(   p.epsz(:) .*exudZoo(:));
exudToON  = sum((1-p.epsz(:)).*exudZoo(:));

mortBac = (fT_h) .* p.mb  .* bac;
mortPhy = (fT_h) .* p.mp  .* phy;
mortZoo = (fT_h) .* (p.mz1 .* max(0,zoo-0.01) + p.mz2 .* max(0,zoo-0.01).^2);

C = 0;
if t > p.lysisTime && p.lysisTime > 0
      C = (1/2 + (1/2)*sin(-pi/2 + (min(1,(iTime-p.lysisTime)/5)*pi)));
      mortPhy(Jehux) = mortPhy(Jehux) + (fT_h) .* 0.5 * C .* phy(Jehux) .* phy(Jehux);
      mortPhy(Jnano) = mortPhy(Jnano) + (fT_h) .* 0.1 * C .* phy(Jnano) .* phy(Jnano);
end

mortTot   = sum(mortBac(:)) + sum(mortPhy(:)) + sum(mortZoo(:));
mortToDIN =      p.omeM  * mortTot;
mortToON  = (1 - p.omeM) * mortTot;

% Infected cells produce carbon compounds that promote bacterial growth.
% These compounds are not explicitly represented in our model.
leak_fac = 1.0;
%D = 0;
if t > (p.lysisTime-2.5) && p.lysisTime > 0
    D = (1/2 + (1/2)*sin(-pi/2 + (min(2,(iTime-p.lysisTime+2.5)/5)*pi)));
    ON_enhanced = ON * (1 + D*Ehux/6);
    lim_ON = ON_enhanced / (ON_enhanced + p.Kon);
    % Extra DMSP and DMS leakage 
    leak_fac = (1 + D*p.leak_inf_fac); % Only in scenario S4
end

% ON uptake by bacteria
uptake_ON_bac  = (UV_bac*fT_h) * p.mub0 * bac * lim_ON;
ON_to_bac      =    p.omeB  * sum(uptake_ON_bac);
remin_ON       = (1-p.omeB) * sum(uptake_ON_bac);

% Loss term to close the balance
sink_ON = p.wON * ON;

dpicodt = uptake_DIN_phy(Jpico) - mortPhy(Jpico) - grazZoo(Jpico+1);
dnanodt = uptake_DIN_phy(Jnano) - mortPhy(Jnano) - grazZoo(Jnano+1);
dEhuxdt = uptake_DIN_phy(Jehux) - mortPhy(Jehux) - grazZoo(Jehux+1);
dzoodt  = preyToZoo - mortZoo;
dbacdt  = ON_to_bac - mortBac - grazZoo(1);

dDINdt  = exudToDIN + mortToDIN + remin_ON - uptake_DIN_phy_tot + input_DIN(min(t,24));
dONdt   = exudToON  + mortToON  - uptake_ON_bac - sink_ON       + input_ON(min(t,24));

%..........................................................................
% DMS(P) model, dependent on the biomass model
%..........................................................................

% DMS loss by photolysis and emission
if p.inputCode == 0
    DMS_sink = p.DMS_loss0 * DMS;
elseif p.inputCode > 0
    DMS_sink = p.DMS_loss1 * DMS;
end

% DMSP:N ratio of phytoplankton cells (vector)
DMSP_ratio = p.thetaSN_P .* (1-(1-p.stress_fac).*lim_DIN);

DMSPp = sum(DMSP_ratio.*1000.*phy);

% Source of DMS and DMSPd at grazing and mortality
graz_to_DMSP  = p.graz_DMSP_frac .* (1000*grazZoo(2:end)) .* DMSP_ratio;
graz_to_DMS   = p.graz_DMS_frac  .* (1000*grazZoo(2:end)) .* DMSP_ratio;
mort_to_DMSP  = p.mort_DMSP_frac .* (1000*mortPhy)        .* DMSP_ratio;
mort_to_DMS   = p.mort_DMS_frac  .* (1000*mortPhy)        .* DMSP_ratio;

% DMSP exudation and DMS leakage (Fraction of nutrient uptake)
leak_DMSP = p.leak_DMSP_frac .* DMSP_ratio .* (1000*uptake_DIN_phy);
leak_DMS  = p.leak_DMS_frac  .* DMSP_ratio .* (1000*uptake_DIN_phy);

% Extra DMSP and DMS release during viral infection
leak_DMSP(4) = leak_DMSP(4) .* leak_fac;
leak_DMS(4)  = leak_DMS(4)  .* leak_fac;

leak_DMSP = leak_DMSP * UV_fac;
leak_DMS = leak_DMS * UV_fac;

% DMPS and DMS limitations
lim_DMS   = DMS   / (p.KDMS   + DMS  );
lim_DMSPd = DMSPd / (p.KDMSPd + DMSPd);

% DMS uptake, depending on model options
Udms_coMet = p.chi_coMet * p.muDMS0_coMet * lim_DMS * (1000*bac);
FDC_diff = 0;
Udms_spec = 0;
% In scenario S4, the proportion of specialist DMS-consumers vary over time
if strcmp(keyVarDMSUptake,'yes')
    alphaB_eff = exp(FDC) / (1+exp(FDC)); 
    Udms_spec = alphaB_eff * p.muDMS0_spec * lim_DMS * (1000*bac);
    FDC_diff = p.nu0 * lim_ON * (1+0.6*C) * (lim_DMS - 0.5 - alphaB_eff*10);
end
Udms = Udms_coMet + Udms_spec;

if strcmp(keyGrowthDMS,'yes')
    Udms = Udms * lim_ON;
end

% Dissolved DMSP uptake, depending on model options
Udmsp = p.muDMSPd0 * lim_DMSPd * (1000*bac);
if strcmp(keyGrowthDMSP,'yes')
    Udmsp = Udmsp * lim_ON;
end
bac_to_DMS = p.alpha2 * Udmsp;

% DMSP uptake by picophytoplankton
Udmsp_phy = lim_DMSPd .* p.muDMSPd0_phy(:) .* p.mup0 .* phy * 1000;
              
dDMSdt   = sum(graz_to_DMS(:))  + sum(mort_to_DMS(:))  + bac_to_DMS + sum(leak_DMS(:)) ...
         - Udms  - DMS_sink;
dDMSPddt = sum(graz_to_DMSP(:)) + sum(mort_to_DMSP(:))              + sum(leak_DMSP(:)) ...
         - Udmsp - sum(Udmsp_phy(:));
dFDCdt    = FDC_diff;

%..........................................................................
% Trends
%..........................................................................

Vdot = [dpicodt;dnanodt;dEhuxdt;dzoodt;dbacdt;dDINdt;dONdt;dDMSdt;dDMSPddt;dFDCdt];

sink_ON = sum(sink_ON);

% Write output fluxes
if strcmp(p.outputFluxes,'yes')
    DMOS_BLOOM_write_fluxes(icounter,p.isim,uptake_DIN_phy,mortPhy,ON_to_bac,remin_ON,mortBac,grazZoo,preyToZoo,mortZoo,exudZoo,sink_ON,...
        leak_DMS,bac_to_DMS,mort_to_DMS,graz_to_DMS,leak_DMSP,mort_to_DMSP,graz_to_DMSP,DMS_sink,Udms,Udmsp,Udmsp_phy) % DMSPp
end

return
