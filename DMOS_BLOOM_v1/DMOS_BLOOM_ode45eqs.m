function [Vdot] = DMOS_BLOOM_ode45eqs(iTime,V0,varargin)
global icounter 
global deltat
global Jpico Jnano Jehux Jzoo Jbac Jdin Jon 
global Jdmspd Jdms Jfdc
global jday
global keyGrowthDMSP keyGrowthDMS keyVarDMSUptake
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
input_DIN = zeros(24,1);
% Add ON input in the fjord
input_ON  = zeros(24,1);
if p.inputCode == 0
    input_ON(:,end)  = sum(p.ON0) * p.wON;
elseif p.inputCode == 1
    input_DIN(1:8)   = input_DIN(1:8) + 1.6;
    input_DIN(14:18) = input_DIN(14:18) + 1.6;
end

%%%%%%%%%%%%%%%%%%%
%DAY OF SIMULATION:
%%%%%%%%%%%%%%%%%%%
%[jday,newday] = DMOS_BLOOM_daycounter(jday,iTime);
jday = DMOS_BLOOM_daycounter(jday,iTime);
t = ceil(iTime);

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
% Plankton model (Le Gland, 29/11/2021)
%..........................................................................

phy = [pico;nano;Ehux];

uptake_DIN_phy = p.mup0 .* phy .* DIN ./ (DIN + p.Kn);
uptake_DIN_phy_tot = sum(uptake_DIN_phy(:));

lim_DIN = DIN ./ (DIN + p.Kn );
lim_ON  = ON  ./ (ON  + p.Kon);

prey = [bac;phy];
% One zooplankton for each prey
grazZoo   = p.muz0 .* zoo .* prey ./ (p.Kz + prey); % Holling type II
exudZoo   = (1 - p.betaz) .* grazZoo;
preyToZoo =      p.betaz   .*grazZoo;  
exudToDIN = sum(   p.epsz(:) .*exudZoo(:));
exudToON  = sum((1-p.epsz(:)).*exudZoo(:));

mortBac = p.mb  .* bac;
mortPhy = p.mp  .* phy;
mortZoo = p.mz1 .* max(0,zoo-0.01) + p.mz2 .* max(0,zoo-0.01).^2;

C = 0;
if t > p.lysisTime && p.lysisTime > 0
      C = (1/2 + (1/2)*sin(-pi/2 + (min(1,(iTime-p.lysisTime)/5)*pi)));
      mortPhy(Jehux) = mortPhy(Jehux) + 0.5 * C .* phy(Jehux) .* phy(Jehux);
      mortPhy(Jnano) = mortPhy(Jnano) + 0.1 * C .* phy(Jnano) .* phy(Jnano);
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
uptake_ON_bac  = p.mub0 * bac * lim_ON;
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
elseif p.inputCode == 1
    DMS_sink = p.DMS_loss1 * DMS;
end

% Amount of DMSP in phytoplankton cells
DMSP_ratio = p.thetaSN_P .* (1-(1-p.stress_fac).*lim_DIN); % vector

%DMSP_ratio(4) = DMSP_ratio(4) * (1 + 2*C); 

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

% DMPS and DMS limitations
lim_DMS   = DMS   / (p.KDMS   + DMS  );
lim_DMSPd = DMSPd / (p.KDMSPd + DMSPd);

% DMS uptake, depending on model options (Le Gland, 18/02/2022)
Udms_coMet = p.chi_coMet * p.muDMS0_coMet * lim_DMS * (1000*bac);
FDC_diff = 0;
Udms_spec = 0;
if strcmp(keyVarDMSUptake,'yes')
    alphaB_eff = exp(FDC) / (1+exp(FDC)); % 
    Udms_spec = alphaB_eff * p.muDMS0_spec * lim_DMS * (1000*bac);
    FDC_diff = p.nu0 * lim_ON * (1+0.5*C) * (lim_DMS - 0.5 - alphaB_eff*10);
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

% Output fluxes for the final round
if strcmp(p.outputFluxes,'yes')
    DMOS_BLOOM_write_fluxes(icounter,p.isim,uptake_DIN_phy,mortPhy,ON_to_bac,remin_ON,mortBac,grazZoo,preyToZoo,mortZoo,exudZoo,sink_ON,...
        leak_DMS,bac_to_DMS,mort_to_DMS,graz_to_DMS,leak_DMSP,mort_to_DMSP,graz_to_DMSP,DMS_sink,Udms,Udmsp,Udmsp_phy) % DMSPp
end

return
