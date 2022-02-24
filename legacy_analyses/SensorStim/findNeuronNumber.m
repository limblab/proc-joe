function [ nn ] = findNeuronNumber(cds, channel, ID, elec )
% returns the neuron number(nn) corresponding to channel, id, elec in cds.
% if it does not exist, returns -1

nn = -1;

i = 1;
while i <= length(cds.units)
    if(cds.units(i).chan == channel && cds.units(i).ID == ID && strcmp(cds.units(i).label, elec))
        nn = i;
        i = 1000;
    end
    i=i+1;
end

end

