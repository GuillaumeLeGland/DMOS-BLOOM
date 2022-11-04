function [day,Ehux,nano,pico,Syn,bac,EhV,MpV,cil,thr,IEC,NO3,Chl,DMS,DMSPd_raw,DMSPd,DMSPp] = DMOS_BLOOM_load_data()

%========================================================================

ncol  = 23; % Number of columns (variables in the data set)
nvar  = 16; % Number of variables
nbag  = 8;  % Number of experiments (including the fjord)
ndat  = 50; % Maximum number of measurements

% Load the external forcings from the input files
fid = fopen('DATA/BERGEN-DATA.csv','r');
fgetl(fid);
D=textscan(fid,repmat('%s ',1,ncol), 'delimiter', ';');
fclose(fid);

index_var = [2,3,4,5,6,7,8,9,11,12,13,14,15,16,18,21]; % List of columns to be read

data = DMOS_BLOOM_read(D,nvar,nbag,ndat,index_var);

day   = squeeze(data(1,:,:));
Ehux  = squeeze(data(2,:,:));
nano  = squeeze(data(3,:,:));
pico  = squeeze(data(4,:,:));
Syn   = squeeze(data(5,:,:));
bac   = squeeze(data(6,:,:));
EhV   = squeeze(data(7,:,:));
MpV   = squeeze(data(8,:,:));
cil   = squeeze(data(9,:,:));
thr   = squeeze(data(10,:,:));
IEC   = squeeze(data(11,:,:));
NO3   = squeeze(data(12,:,:));
Chl   = squeeze(data(13,:,:));
DMS   = squeeze(data(14,:,:));
DMSPd_raw = squeeze(data(15,:,:));
DMSPp = squeeze(data(16,:,:));

% Convert cels/mL to uM-N
Ehux = Ehux * 14.4 / (12000*6.625); % 14.4 pg of carbon per cell, following Vincent et al. (2021) 
nano = nano * 14.4 / (12000*6.625);
pico = pico * 1.8  / (12000*6.625); % 1.8 pgC per cell, ased on a sphere with a radius of 1.25µm (Micromonas)
Syn  = Syn  * 0.1  / (12000*6.625); % 0.1 pgC per cell
bac  = bac  * 0.01 / (12000*5.1); % 0.01 pgC per cell, following the supplement of Vincent et al. [2021]

cil  = cil  * 1000 / (12000*5.5); % 1000 pgC per cell (based on the volume of 5000 µm3 observed by flow cameras)
% On average, thrautochytrids have 165 pg of carbon per cell (Kimura et al., 1999) and 181 18s copies per cell
thr  = thr  * 165  / (12000*5.1*181);

% mcp copies / mL -> mcp copies / L
EhV   = EhV * 1000;

% Eliminate outliers that are not shown in Vincent et al 2022
bac(1,11)  = NaN;
bac(1,22) = NaN;
bac(3,26) = NaN;
bac(4,31) = NaN;

% Smoothed version of dissolved DMSP
% First remove sawtooth spikes (>40nM with all neighbors <40nM) 
% and replace them with interpolates.
% Then cap all DMSPd values by 20% of total DMSP. This will remove fjord 
% outliers contaminated by raft algae fromday 15 onward.
% Care! at day 1 and 17, DMSPd may really be between 50 and 80. This
% result is consistent between bags and do not cross the 20% threshold
% because particulate DMSP is high at that time
num_bag = size(Ehux,1);
tmax = size(DMSPd_raw,2);
DMSPd_smooth = DMSPd_raw;
DMSPp_smooth = DMSPp;
for b=1:num_bag
    id = ~isnan(DMSPd_raw(b,:));
    idmax = sum(id(:));
    dat = zeros(idmax,1);
    ji = 1;
    for t=1:tmax
        if id(t) == 1
            dat(ji) = t;
            ji = ji + 1;
        end
    end
    for ji=2:idmax-1;
        if ((DMSPd_raw(b,dat(ji)) > 40 && DMSPd_raw(b,dat(ji-1)) < 40 && DMSPd_raw(b,dat(ji+1)) < 40) ...
                || DMSPd_raw(b,dat(ji)) > 100) ...
                || (b == 8 && (dat(ji) == 33 || dat(ji) == 37) )
            c1 = (day(b,dat(ji+1))-day(b,dat(ji))  )/(day(b,dat(ji+1))-day(b,dat(ji-1)));
            c2 = (day(b,dat(ji))  -day(b,dat(ji-1)))/(day(b,dat(ji+1))-day(b,dat(ji-1)));
            swap = DMSPd_raw(b,t);
            DMSPd_smooth(b,dat(ji)) = c1*DMSPd_raw(b,dat(ji-1)) + c2*DMSPd_raw(b,dat(ji+1));
            DMSPp_smooth(b,t) = DMSPp(b,t) + swap - DMSPd_smooth(b,t); 
        end
    end
    if DMSPd_raw(b,dat(idmax)) > 40
       tarr = 1:tmax;
       d = max(tarr(DMSPd_raw(b,:)<40));
       e = DMSPd_raw(b,d);
       for t=d+1:tmax
           if DMSPd_raw(b,t) > 40
               DMSPd_smooth(b,t) = e;
           end
       end
    end
end
for b=1:num_bag
    for t=1:tmax
        if DMSPd_smooth(b,t) > 0.2*(DMSPd_smooth(b,t)+DMSPp(b,t))
            swap = DMSPd_smooth(b,t); 
            DMSPd_smooth(b,t) = 0.2*(DMSPd_smooth(b,t)+DMSPp(b,t));
            DMSPp_smooth(b,t) = DMSPp_smooth(b,t) + swap - DMSPd_smooth(b,t);
        end
    end
end
% Bag 6, day 16 is likely goood, by comparison with other bags
DMSPd_smooth(6,33) = DMSPd_raw(6,33);
DMSPp_smooth(6,33) = DMSPp(6,33);
% In fjord, data after day 22 are likely contaminated
DMSPd_smooth(8,47:49) = DMSPd_smooth(8,46);
DMSPp_smooth(8,47:49) = DMSPp_smooth(8,46);
DMSPd = DMSPd_smooth;
DMSPp = DMSPp_smooth;

return

function [data_in] = DMOS_BLOOM_read(D,nvar,nbag,ndat,index_var)

data_in  = zeros(nvar,nbag,ndat);

b = 1;
j = 0;
for i = 1:length(D{1})
    bn = str2double(char(D{1}(i)));% Bag number (new)
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

