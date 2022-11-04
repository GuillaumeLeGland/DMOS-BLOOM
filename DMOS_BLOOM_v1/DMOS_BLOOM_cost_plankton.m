function [ res ] = DMOS_BLOOM_cost_plankton(iTime,X,X0,Xtol,p,pOpt,day,obs,obs_tol)
%DMOS_BLOOM_cost_plankton
%   Cost function of the SUMMIT_BERGEN model for plankton 
%   Inputs are the unknown parameters to be optimized
%   Ouputs are the differences between model and observed plankton
%   concentrations

%--------------------------------------------------------------------------
% Simulations (Fjord + standard bag)
%--------------------------------------------------------------------------

% Not effective if two parameters share the same value -> solved
% Still not effective if combination of various traits -> can be solved in
% the model itself
% paramopt can therefore be larger than X
nx = length(X);
fnVar = fieldnames(pOpt);
nVar = numel(fnVar);
for ji=1:nVar
    arr = pOpt.(fnVar{ji});
    sz = length(arr);
    for n=1:sz
        if arr(n) > 0
            p.(fnVar{ji})(n) = X(arr(n));
        end
    end
end
%p.Kon = p.mub0 / p.Kon;
%p.Kz  = p.muz0 ./ p.Kz;

cost_parms = nanmean(((X-X0)./Xtol).^2);
res_parms = (X - X0) ./ (Xtol * sqrt(nx));
cost_parms_bis = nansum(res_parms.^2);

V0 = [p.pico0;
      p.nano0;
      p.Ehux0;
      p.zoo0;
      p.bac0;
      p.DIN0;
      p.ON0;
      p.dmspd0;
      p.dms0;
      p.fdc0];

% Remove fields that the model does not need
%p = rmfield(p,{'pico0','nano0','Ehux0','zoo0','bac0',...
%    'DIN0','ON0','dmspd0','dms0','bs0'});
p = rmfield(p,{'pico0','nano0','Ehux0','zoo0','bac0',...
    'DIN0','dmspd0','dms0','fdc0'});

%dbstop;

nsteps = length(iTime);
deltat = iTime(2) - iTime(1);
nout = length(V0);
Vout = zeros(nsteps,nout,2);

% 1) Control without DIN input -> Fjord
p.inputCode = 0;
p.lysisTime = -1;
%p.lysisTime = 20;
Vout(:,:,1) = ode4(@DMOS_BLOOM_ode45eqs,iTime,V0,p);

% 2) Standard with DIN input -> Bags 1, 2, 3, 5 and 6
p.inputCode = 1;
%p.lysisTime = 20; % Viral mortality at the end,even in bags 1-2-3-5-6 (Le Gland, 07/03/2022) 
Vout(:,:,2) = ode4(@DMOS_BLOOM_ode45eqs,iTime,V0,p);

%--------------------------------------------------------------------------
% Extraction of model results at points where observations exist
%--------------------------------------------------------------------------

%id_var = [1,2,3,5,6,7]; % pico, nano, Ehux, bac, cil and DIN are compared to observations
id_var = [1,2,4,9,10,11]; % pico, nano, Ehux, Ehux pred, bacteria and DIN are compared to observations

mod = zeros(size(obs));

for ji = 1:size(obs,2)
    d = day(ji);
    t = d/deltat + 1;
    for jj = 1:size(obs,3)
        for js = 1:size(obs,1);
            if isnan(obs(js,ji,jj)) || (jj == 6 && d > 14)
                mod(js,ji,jj) = NaN; %No data
            else
                id = id_var(jj);
                if id == 2
                    nanoconc = Vout(min(t,end),2:3,js); % Sum of the two nanophytoplanktons
                    mod(js,ji,jj) = sum(nanoconc(:)); 
                else
                    mod(js,ji,jj) = Vout(min(t,end),id,js);
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
% Cost function and residuals
%--------------------------------------------------------------------------

% Cost of distance between model and observations
cost_dist = nanmean( ( (mod(:) - obs(:)) ./ obs_tol(:) ) .^ 2);
num_dist  = sum(mod(:) == mod(:));
res_dist  = (mod(:) - obs(:)) ./ (obs_tol(:)*sqrt(num_dist));
cost_dist_bis = nansum(res_dist.^2);

%cost = cost_parms + cost_ss + cost_dist;
cost = cost_parms + cost_dist;

res = [res_parms;res_dist];
%res = [res_parms;res_ss;res_dist];
res(res~=res) = 0.0; % Zero out NaN values

cost_bis = sum(res.^2);

%dbstop;

disp('-------------------- DMOS_BLOOM_cost_plankton --------------------')
disp(['V0:  ', mat2str(V0',3)])
disp(['X:  ', mat2str(X',3)])

disp(['mup0:  ', num2str(p.mup0'), '  Kn:  ', num2str(p.Kn')])
disp(['mp:  ', num2str(p.mp')])
disp(['mub0:  ', num2str(p.mub0), '  Kon:  ', num2str(p.Kon'), '  omeB:  ', num2str(p.omeB)])
disp(['mb:  ', num2str(p.mb)])
disp(['muz0:  ', num2str(p.muz0'), '  Kz:  ', num2str(p.Kz')])
disp(['betaz:  ', num2str(p.betaz'), '  epsz:  ', num2str(p.epsz')])
disp(['mz1:  ', num2str(p.mz1'), '  mz2:  ', num2str(p.mz2')])
disp(['omeM:   ', num2str(p.omeM), '  wON:  ', num2str(p.wON)])
  
disp(['cost:    ', num2str(cost)])
disp(['cost_parms:   ', num2str(cost_parms), '   cost_dist:   ', num2str(cost_dist)]) %, ...
 %'   cost_ss:   ', num2str(cost_ss)])
disp(['cost_bis:    ', num2str(cost_bis)])
disp(['cost_parms_bis:   ', num2str(cost_parms_bis), '   cost_dist_bis:   ', num2str(cost_dist_bis)]) %, ...
 %'   cost_ss_bis:   ', num2str(cost_ss_bis)])
disp('---------------------------------------------------------------------')
  
%dbstop

return