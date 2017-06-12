function [ prob ] = probSpike( current, distances, threshold, gainThreshold, s )
% computes the probability of an ICMS pulse evoking a pulse
stdThreshold = s*threshold;
current = current./distances.^2;
Iz = (gainThreshold*current - threshold)./stdThreshold;
prob = normcdf(Iz,0,1); % ~ N(0,1)

end

