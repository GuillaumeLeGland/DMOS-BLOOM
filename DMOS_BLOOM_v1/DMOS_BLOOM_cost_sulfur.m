function [ res ] = DMOS_BLOOM_cost_sulfur(iTime,Y,Y0,Ytol,p,pOpt,day,obs,obs_tol)
%SUMMIT_BERGEN_cost
%   Cost function of the DMOS_BLOOM model for DMS(P)
%   Inputs are the unknown parameters to be optimized
%   Ouputs are the differences between model and observed plankton
%   concentrations

%--------------------------------------------------------------------------
% Simulations (Fjord + standard bag)
%--------------------------------------------------------------------------

ny = length(Y);
fnVar = fieldnames(pOpt);
nVar = numel(fnVar);
for ji=1:nVar
    arr = pOpt.(fnVar{ji});
    sz = length(arr);
    for n=1:sz
        if arr(n) > 0
            p.(fnVar{ji})(n) = Y(arr(n));
        end
    end
end
%dbstop;

% Considerations of mass balance are introduced here
p.mort_DMSP_frac = 1 - p.mort_DMS_frac;
graz_frac = 0.7;
p.graz_DMS_frac = graz_frac * p.mort_DMS_frac;
p.graz_DMSP_frac = graz_frac - p.graz_DMS_frac;

cost_parms = mean(((Y-Y0)./Ytol).^2);
res_parms = (Y - Y0) ./ (Ytol * sqrt(ny));
cost_parms_bis = nansum(res_parms.^2);

V0 = [p.pico0;
      p.nano0;
      p.Ehux0;
      p.zoo0;
      p.bac0;
      p.DIN0;
      p.ON0;
      p.dms0;
      p.dmspd0;
      p.fdc0];

% Remove fields that the model does not need
p = rmfield(p,{'pico0','nano0','Ehux0','zoo0','bac0',...
    'DIN0','dmspd0','dms0','fdc0'});

nsteps = length(iTime);
deltat = iTime(2) - iTime(1);
nout =length(V0);
Vout = zeros(nsteps,nout,2);

% 1) Control without DIN input -> Fjord
p.inputCode = 0;
p.lysisTime = -1;
Vout(:,:,1) = ode4(@DMOS_BLOOM_ode45eqs,iTime,V0,p);

% 2) Standard with DIN input -> Bags 1, 2, 3, 5 and 6
p.inputCode = 1;
Vout(:,:,2) = ode4(@DMOS_BLOOM_ode45eqs,iTime,V0,p);

%--------------------------------------------------------------------------
% Extraction of model results at points where observations exist
%--------------------------------------------------------------------------

id_var = [13,14]; % DMS and DMSPd are compared to observations
obs = obs(:,:,1:2);
obs_tol = obs_tol(:,:,1:2);
mod = zeros(size(obs));

for ji = 1:size(obs,2)
    d = day(ji);
    t = d/deltat + 1;
    for jj = 1:2 %1:size(obs,3)
        for js = 1:size(obs,1);
            if isnan(obs(js,ji,jj)) % || (jj == 6 && d > 14)
                mod(js,ji,jj) = NaN; %No data
            else
                id = id_var(jj);
                mod(js,ji,jj) = Vout(min(t,end),id,js);
            end
        end
    end
end

%--------------------------------------------------------------------------
% Cost function and residuals
%--------------------------------------------------------------------------

% Cost of distance between model and observations
cost_dist = nanmean( ( (mod(:) - obs(:)) ./ obs_tol(:) ) .^ 2);
num_dist  = sum(obs(:) == obs(:));
res_dist  = (mod(:) - obs(:)) ./ (obs_tol(:)*sqrt(num_dist));
cost_dist_bis = nansum(res_dist.^2);

cost = cost_parms + cost_dist;

res = [res_parms;res_dist];
res(res~=res) = 0.0; % Zero out NaN values

cost_bis = sum(res.^2);

disp('-------------------- DMOS_BLOOM_cost_sulfur --------------------')
disp(['V0:  ', mat2str(V0',3)])
disp(['Y:  ', mat2str(Y',3)])

disp(['DMS_loss:  ', num2str([p.DMS_loss0, p.DMS_loss1]), '  muDMS0_coMet:  ', num2str(p.muDMS0_coMet), '  KDMS:  ', num2str(p.KDMS)])
disp(['muDMS0_spec:  ', num2str(p.muDMS0_spec), '  nu0:  ', num2str(p.nu0)])
disp(['muDMSPd0:  ', num2str(p.muDMSPd0), '  KDMSPd:  ', num2str(p.KDMSPd)])
disp(['muDMSPd0_phy:  ', num2str(p.muDMSPd0_phy')])
disp(['alpha2:  ', num2str(p.alpha2), '  chi_coMet:  ', num2str(p.chi_coMet)])
disp(['stress_fac:  ', num2str(p.stress_fac), '  thetaSC_P:  ', num2str(p.thetaSN_P')])
disp(['graz_DMSP_frac:  ', num2str(p.graz_DMSP_frac'), '  graz_DMS_frac:  ', num2str(p.graz_DMS_frac')])
disp(['mort_DMSP_frac:  ', num2str(p.mort_DMSP_frac'), '  mort_DMS_frac:  ', num2str(p.mort_DMS_frac')])
disp(['exud_DMSP_frac:  ', num2str(p.leak_DMSP_frac'), '  exud_DMS_frac:  ', num2str(p.leak_DMS_frac')])
  
disp(['cost:    ', num2str(cost)])
disp(['cost_parms:   ', num2str(cost_parms), '   cost_dist:   ', num2str(cost_dist)])
disp(['cost_bis:    ', num2str(cost_bis)])
disp(['cost_parms_bis:   ', num2str(cost_parms_bis), '   cost_dist_bis:   ', num2str(cost_dist_bis)])
disp('---------------------------------------------------------------------')
  
%dbstop

return