function [day_Meth,Meth_free,Meth_part] = DMOS_BLOOM_load_DMS_consumers()
%DMOS_BLOOM_load_DMS_consumers
%   Load relative abundances of Methylophaga from input file
%   Other DMS consumers are also present in the data

%========================================================================

ncol = 14;

% Load the external forcings from the input files
fid = fopen('DATA/BERGEN-DMS-CONSUMERS.csv','r');
fgetl(fid);
D=textscan(fid,repmat('%s ',1,ncol), 'delimiter', ';');
fclose(fid);

nvar  = 3; % Number of variables
nbag  = 8;  % Number of experiments (including the fjord)
ndat  = 25; % Maximum number of measurements

index_var = [2,3,8]; % List of columns to be read

data = DMOS_BLOOM_read(D,nvar,nbag,ndat,index_var);

day_Meth   = squeeze(data(1,:,:));
Meth_free  = squeeze(data(2,:,:));
Meth_part  = squeeze(data(3,:,:));

return

function [data_in] = DMOS_BLOOM_read(D,nvar,nbag,ndat,index_var)

data_in  = zeros(nvar,nbag,ndat);

b = 1;
j = 0;
for i = 1:length(D{1})
    bn = str2double(char(D{1}(i))); % Bag number (new)
    if bn > b
        j = 1;
        b = bn;
    else
        j = j + 1;
    end
    for v = 1:nvar
        id = index_var(v);
        data_in(v,b,j) = str2double(char(D{id}(i)));
    end
end

return