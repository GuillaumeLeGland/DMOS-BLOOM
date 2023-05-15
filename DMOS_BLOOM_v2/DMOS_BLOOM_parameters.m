function [param,paramOptX,paramOptY,X0,Y0,xtl,ytl,Xtol,Ytol,lowBndX,uppBndX,lowBndY,uppBndY] = DMOS_BLOOM_parameters()

%========================================================================
% Define all numerical, physical and biological constants of the model
%========================================================================

%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL CONSTANTS
%%%%%%%%%%%%%%%%%%%%%
%========================================================================
%........................................................................

% Select one of the four scenarios described in the article (1/2/3/4)
% You can also change each parameter by hand
% Scenario = 1;
% Scenario = 2;
% Scenario = 3;
Scenario = 4;

param.outputFluxes = 'not';

param.ndays = 24; % Number of days in simulation
param.deltat  = 1/24; % Time step (d)

%..........................................................................
% Microbial model parameters
%..........................................................................

% Physical parameters
param.temp = [11.52,11.52,11.52,11.58,12.67,12.88,13.33,13.81,14.01,15.68,16.28,16.37,14.49,13.22,14.17,15.12,15.82,...
              16.48,15.04,12.91,11.80,11.27,10.88,11.27;...
              11.48,11.99,11.91,11.88,12.63,13.08,13.44,14.11,14.78,15.95,16.83,16.81,14.97,14.14,14.66,15.19,15.84,...
              16.51,15.40,13.40,12.20,11.39,10.91,11.34;...
              11.64,11.99,11.87,11.79,12.61,13.06,13.43,14.08,14.74,15.94,16.78,16.69,14.78,14.07,14.61,15.15,15.76,...
              16.45,15.30,13.33,12.06,11.38,10.89,11.30;
              11.28,11.99,11.94,11.71,12.63,13.11,13.48,14.15,14.51,15.99,16.81,16.91,15.13,14.18,14.76,15.35,15.93,...
              16.54,15.61,13.51,12.32,11.38,10.90,11.36];
param.PAR  = [133.11,131.73,135.45,129.35,127.60,131.21,126.82,127.87,129.80,130.02,130.85,139.11,131.72,134.16,136.6,...
    139.04,114.87,59.49,57.13,88.71,111.77,10.28,65.37,115.81];


% PHYTOPLANKTON
param.mup0 = [0.7; 0.7; 0.7; 0.7]; % Maximum growth rate (d-1)
param.Kn   = [0.5; 0.5; 0.5; 0.5]; % Half-saturation constant for DIN (然)
param.mp   = [0.05; 0.05; 0.05; 0.05]; % Linear mortality rate (d-1)

% HETEROTROPHIC BACTERIA
param.mub0 = 5.0;  % Maximum uptake rate (d-1)
param.Kon  = 5.0;  % Half-saturation constant for ON (然)
%param.Aon  = param.mub0 / param.Kon; % Affinity for ON (d-1 然-1)
param.omeB = 0.25; % Assimilation efficiency (the rest goes to DIN) (-)
param.mb   = 0.1;  % Linear mortality rate (d-1)

% ZOOPLANKTON (one for each prey type)
param.muz0  = 1.4 * ones(5,1);  % Maximum grazing rate (d-1)
param.Kz    = 4.0 * ones(5,1);  % Half-saturation constant for food (然-N)
param.Kz(1) = 2.0;    
param.betaz = 0.7 * ones(5,1);  % Assimilation efficiency (-)
param.epsz  = 0.25 * ones(5,1); % Fraction of unassimilated food going to DIN (the rest goes to ON) (-)
param.mz1   = 0.05 * ones(5,1); % Linear mortality rate (d-1)
param.mz2   = 0.5 * ones(5,1);  % Quadratic mortality rate (然-1 d-1)

param.omeM = 0.25; % Fraction of mortality (all plankton types) going to DIN (the rest goes to ON) (-)

%param.wON = 0.07;  % ON removal rate (d-1) due to sinking
param.wON = 0.05; % New value after revision to better reproduce PON observations

%..........................................................................

% Initial nitrate concentration (然)
param.DIN0  = 0.27;
param.ON0   = 1.5;  % Around 10.0 in observations. We will suppose 15% is labile.

% Initial concentration of phytoplankton species, in 然-N (Le Gland, 29/11/2021)
param.pico0 = 0.45;
param.nano0 = [0.3; 0.005];
param.Ehux0 = 0.01;
param.bac0  = 0.1;
param.zoo0 = [0.05;0.225;0.15;0.0025;0.005];

%..........................................................................
% Sulfur model parameters
%..........................................................................

% DMSP:C ratio of phytoplankton
if Scenario >= 2
    thetaSC_P = [0.018; 0.018; 0.011; 0.011]; % Scenarios S2, S3 and S4
    param.stress_fac = 0.2; % DMSP:C is reduced to 20% of maximum under nutrient-replete conditions
else
    thetaSC_P = [0.009; 0.009; 0.009; 0.009]; % Scenario S1
    param.stress_fac = 1.0;
end
param.thetaSN_P = thetaSC_P * 6.625; 

param.thetaSC_B = 0.012;

% Abiotic DMS sink (photolysis + gas exchange)
param.DMS_loss0 = 0.125; % DMS abiotic removal rate (d-1) in the fjord
param.DMS_loss1 = 0.2;   % DMS abiotic removal rate in the mesocosms

%Biotic DMS sink (bacterial uptake)
param.muDMS0_coMet = 0.15; %0.2; % Maximum uptake rate by co-metabolizers (?)
param.muDMS0_spec  = 4.5; %6.0;  % Maximum uptake rate by specialists (?)
param.KDMS         = 8.0; % Half-saturation constant (nM)

% Dissolved DMSP uptake
if Scenario >= 3
    param.muDMSPd0 = 0.55; %0.6;  % Maximum DMSPd uptake rate by bacteria (?)
else
    param.muDMSPd0 = 0.8; %0.9;
end
param.KDMSPd   = 20.0; % Half-saturation constant (nM)
% Maximum uptake rate by phytoplankton
if Scenario >= 3
    param.muDMSPd0_phy = [0.055; 0.0; 0.0; 0.0]; %[0.06; 0.0; 0.0; 0.0]; % Scenarios S3 and S4
else
    param.muDMSPd0_phy = zeros(4,1); % Scenarios S1 and S2
end

% Bacterial DMS yield (DMS production / DMSPd uptake)
if Scenario >= 3
    param.alpha2 = 0.28; % Scenarios S3 and S4
else
    param.alpha2 = 0.16; % Scenarios S1 and S2
end
    
% Fraction of DMS consumers
param.chi_coMet = 0.2; % Constant fraction of DMS co-metabolizers
if Scenario >= 4
    param.chi_spec0 = log(0.0005/0.9995); % Initial fraction of DMS specialists (scenario S4)
    param.nu0 = 2.5; % Maximum competitive advantage of DMS consumers under DMS-replete conditions (d-1)
else
    param.chi_spec0 = log(0.002/0.998); % Scenarios S1, S2 and S3
    param.nu0 = 0.0;
end


% Fate of particulate DMSP at mortality and grazing and leakage rates
% Immediate DMSP cleavage at cell death is neglected
param.graz_DMSP_frac = [0.7; 0.7; 0.7; 0.7];
param.graz_DMS_frac  = [0.0; 0.0; 0.0; 0.0];
param.mort_DMSP_frac = [1.0; 1.0; 1.0; 1.0];
param.mort_DMS_frac  = [0.0; 0.0; 0.0; 0.0];
param.leak_DMSP_frac = [0.1; 0.1; 0.1; 0.1];
param.leak_DMS_frac  = [0.0; 0.06; 0.06; 0.2];

% Extra DMSP and DMS leakage during infection in scenario S4
if Scenario >= 4
    param.leak_inf_fac = 2.5; %2.8; % 280% increase
else
    param.leak_inf_fac = 0.0;
end

% Initial concentration of sulfur coumpounds (nM)
param.dmspd0 = 6.0;
param.dms0   = 8.0; %5.0;

% Initial fraction of specialist DMS-consuming bacteria
param.fdc0 = param.chi_spec0;

%..........................................................................
% Initial value for the vectors of unknowns
%..........................................................................

%Initial value of plankton unknowns
X0 = [param.mup0;
      param.mub0;
      param.Kon;
      param.omeB;
      param.muz0;
      param.Kz;
      param.mz2;];
  
% Initial value of sulfur unknowns
Y0 = [param.muDMS0_coMet;
      param.muDMSPd0;
      param.muDMSPd0_phy(1);
      param.alpha2;
      param.leak_DMS_frac(2:2:4);
      param.leak_DMSP_frac(1);
      param.KDMS;
      param.KDMSPd;
      param.nu0;
      param.mort_DMS_frac(4)];
  
% Labels when axis is list of unknown variables (xtl = plankton / ytl = sulfur)
xtl = {'mup01','mup02','mup03','mup04','mub0','Kon','omeB','muz01','muz02','muz03','muz04','muz05',...
       'Kz1','Kz2','Kz3','Kz4','Kz5','mz21','mz22','mz23','mz24','mz25'};
ytl = {'muDMS0','muDMSPd0','muDMSPd0_phy','alpha2','eDMSfrac23','eDMSfrac4',...
       'eDMSPfrac','KDMS','KDMSPd','nu0','DMSmortFrac4'};

% Tolerance for each unknown
Xtol = X0;
Ytol = Y0;
  
% Linking unknowns to the corresponding parameters (plankton)
paramOptX.mup0  = 1:4; % All four mup0 require optimization
paramOptX.mub0  = 5; 
paramOptX.Kon   = 6;
paramOptX.omeB  = 7;
paramOptX.muz0  = [8,9,10,11,12];  
paramOptX.Kz    = [13,14,15,16,17]; 
paramOptX.mz2   = [18,19,20,21,22];

% Increased tolerance on the poorly known mz2. Reduced tolerance on the 
% phytoplankton growth rates and bacterial assimilation efficiency.
Xtol(18:22) = 2*Xtol(18:22);
Xtol(1:4) = 0.5*Xtol(1:4);
Xtol(7) = 0.5*Xtol(7);

% Lower and upper boundaries for the plankton unknowns
nx = length(X0);
lowBndX = zeros(nx,1);
lowBndX(13) = 0.5;
lowBndX(14:17) = 1.0;
uppBndX = [2 2 2 2 15 15 1 10 10 10 10 10 10 10 10 10 10 5 5 5 5 5];

% Linking unknowns to the corresponding parameters (sulfur)
paramOptY.muDMS0_coMet    = 1;
paramOptY.muDMSPd0  = 2;
paramOptY.muDMSPd0_phy = [3,0,0,0];
paramOptY.alpha2    = 4;
paramOptY.leak_DMS_frac = [0,5,5,6];
paramOptY.leak_DMSP_frac = [7,7,7,7];
paramOptY.KDMS      = 8;
paramOptY.KDMSPd    = 9;
paramOptY.nu0    = 10;
paramOptY.mort_DMS_frac = [0,0,0,11];

Ytol(1:3) = 10*Ytol(1:3);
Ytol(8:9) = 10*Ytol(8:9);
Ytol(11)  = 1;

% Lower and upper boundaries for the sulfur unknowns
ny = length(Y0);
lowBndY = zeros(ny,1);
lowBndY(4)   = 0.05;
lowBndY(8:9) = 8.0;
uppBndY = ones(ny,1);
uppBndY(1:3)  = 100.0;
uppBndY(8:10) = 100.0; 

%dbstop;

%........................................................................
%========================================================================
%************************************************************************

return