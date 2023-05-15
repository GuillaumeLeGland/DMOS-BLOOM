function [jday,newday] = DMOS_BLOOM_daycounter(jday,iTime) 

jdaynew = ceil(iTime); % Current day
if jdaynew == jday
    newday = 'not';
else
    newday = 'yes';
    jday = jdaynew;
end

return
    
