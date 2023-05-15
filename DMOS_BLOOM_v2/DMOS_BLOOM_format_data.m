function [obs,obs_tol] = DMOS_BLOOM_format_data(raw_obs,weight_obs,bag_exp,keyObsTolerance)
% Format observations into array usable in cost functions (Le Gland, 15/12/2021)

nexp  = max(bag_exp); % Number of experiments
[~,npoint,nvar] = size(raw_obs); 

obs = zeros(nexp,npoint,nvar);
obs_tol = zeros(nexp,npoint,nvar);

for ji = 1:nexp
    idx = (bag_exp == ji);
    ens = sum(idx); % Ensemble size
    for jj = 1:nvar
        var_idx = raw_obs(idx,:,jj);
        obs(ji,:,jj) = nanmean(var_idx,1);
        if ens > 1
            obs_tol(ji,:,jj) = max(0.05,nanstd(var_idx,1)) * sqrt(ens/(ens-1)) / weight_obs(jj);
        else
            obs_tol(ji,:,jj) = max(0.2,obs(ji,:,jj)) / weight_obs(jj);
        end
        if strcmp(keyObsTolerance,'global')
            obs_tol(ji,:,jj) = nanmean(obs_tol(ji,:,jj)) / 2;
        end
    end
end

end