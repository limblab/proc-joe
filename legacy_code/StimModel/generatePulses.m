function [pulses] = generatePulses(numPulses, magPulse, pulseWidth, timeBetween)
% this function generates a list of pulses ([time, mag, PW]) based on the
% provided data
pulses = zeros(numPulses,3);
for i = 1:numPulses
    if(i == 1)
        pulses(i,1) = 0;
    else
        pulses(i,1) = pulses(i-1,1) + timeBetween;
    end
    pulses(i,2) = magPulse;
    pulses(i,3) = pulseWidth;
end


end