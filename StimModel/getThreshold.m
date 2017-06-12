function [ threshold ] = getThreshold(time, timeAbs, tau, gamma, pulseWidth, chronaxie, rheobase)
% computes the current threshold for producing a spike based on the time,
% and some parameters

if(time <= 0) % no spike has occurred
    threshold = ones(size(timeAbs,1),1)*rheobase*(1+chronaxie/pulseWidth);
elseif(time < timeAbs) % spike has occurred within the absolute refractory period
    threshold = 1000000;
else % relative refractory period
    initThreshold = rheobase*(1+chronaxie/pulseWidth);
    threshold = initThreshold*(1 + gamma*exp(-(time-timeAbs)/tau));
end

end

