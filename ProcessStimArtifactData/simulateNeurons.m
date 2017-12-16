%% this function takes a cds with artifactData, adds artificial neurons to
% channels without units (determined by cds or by user settings), filters
% and thresholds those artifact snippets and generates a NEV file. Then,
% the simulated neurons can be sorted, cds generated, and a reliability
% metric can be determined.

%% settings
chansUse = [];
% chansUse = [28,30,31,32,33,34,36];
maxChansUse = 7;
threshold = -3.5;
idxsAddNeurons = [20:20:150];

%% add artificial neurons to artifact data on specified channels
if(isempty(chansUse)) % find channels based on unit data
    chansUse = 1:96;
    chansWithUnits = [];
    for i = 1:size(cds.units,2)
        if(cds.units(i).ID == 1)
            chansWithUnits(end+1) = cds.units(i).chan;
        end
    end
    chansWithUnits = unique(chansWithUnits);
    chansUse = setdiff(chansUse,chansWithUnits);
    chansUse = sort(datasample(chansUse,maxChansUse,'replace',false));
end

