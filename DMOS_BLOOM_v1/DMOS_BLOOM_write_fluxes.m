function DMOS_BLOOM_write_fluxes(icounter,isim,uptake_DIN_phy,mort_phy,ON_to_bac,remin_ON,mortBac,graz_zoo,preyToZoo,mortZoo,exud,sink_ON,...
    leak_DMS,bac_to_DMS,mort_to_DMS,graz_to_DMS,leak_DMSP,mort_to_DMSP,graz_to_DMSP,DMS_sink,Udms,Udmsp,Udmsp_phy) % DMSPp

%========================================================================

global FLUXES

%--------------------------------------------------------------------------
% Plankton model FLUXESs
%--------------------------------------------------------------------------

FLUXES.uptake_DIN_phy(:,icounter,isim) = uptake_DIN_phy;
FLUXES.mort_phy(:,icounter,isim) = mort_phy;

FLUXES.ON_to_bac(icounter,isim) = ON_to_bac;
FLUXES.remin_ON(icounter,isim) = remin_ON;
FLUXES.mort_bac(icounter,isim) = mortBac;

FLUXES.graz_zoo(:,icounter,isim) = graz_zoo;
FLUXES.prey_to_zoo(:,icounter,isim) = preyToZoo;
FLUXES.mort_zoo(:,icounter,isim) = mortZoo;
FLUXES.exud(:,icounter,isim) = exud;

FLUXES.sink_ON(icounter,isim) = sink_ON;

%--------------------------------------------------------------------------
% Sulfur model FLUXESs
%--------------------------------------------------------------------------

FLUXES.leak_DMS(:,icounter,isim) = leak_DMS;
FLUXES.bac_to_DMS(icounter,isim) = bac_to_DMS;
FLUXES.mort_to_DMS(:,icounter,isim) = mort_to_DMS;
FLUXES.graz_to_DMS(:,icounter,isim) = graz_to_DMS;

FLUXES.leak_DMSP(:,icounter,isim) = leak_DMSP;
FLUXES.mort_to_DMSP(:,icounter,isim) = mort_to_DMSP;
FLUXES.graz_to_DMSP(:,icounter,isim) = graz_to_DMSP;

FLUXES.DMS_sink(icounter,isim) = DMS_sink;
FLUXES.Udms(icounter,isim) = Udms;

FLUXES.Udmsp(icounter,isim) = Udmsp;
FLUXES.Udmsp_phy(:,icounter,isim) = Udmsp_phy;

% FLUXES.DMSPp(icounter,isim) = DMSPp;

return

