function [keyInvPlankton,keyInvSulfur,keyObsTolerance,keyGrowthDMSP,keyGrowthDMS,keyVarDMSUptake] = DMOS_BLOOM_keys()

%========================================================================
% Keys to activate different options
%........................................................................
% Estimate plankton parameters by inversion ('yes') or use prescribed values ('not')
keyInvPlankton = 'not';
% Same for DMS(P) parameters
keyInvSulfur = 'not';
% The weight of each observation in inversion is based on local standard deviation 
% ('local') or on standard deviation averaged over the whole experiment ('global')
keyObsTolerance = 'global';
%........................................................................
% DMSP uptake proportional to growth (yes/not)
keyGrowthDMSP = 'yes';
% DMS  uptake proportional to growth (yes/not)
keyGrowthDMS  = 'yes'; 
% Proportion of DMS-consuming bacteria vary over time (yes/not)
keyVarDMSUptake = 'yes';   
%........................................................................
%========================================================================
%************************************************************************
return


