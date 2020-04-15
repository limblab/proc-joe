%% this function removes the unit ID from a cds object and replaces all sorting with unsorted
units = cds.units;
units_comb = [];
for chan = 1:96
    unit_idx = find([units.chan] == chan);
    spikes = [];
    for u = unit_idx
        spikes = [spikes;units(u).spikes];
    end
    
    if(~isempty(unit_idx))
        [~,sort_idx] = sort(spikes.ts);
        spikes = spikes(sort_idx,:);
        
        if(isempty(units_comb))
            units_comb = units(unit_idx(1));
        else
            units_comb(end+1) = units(unit_idx(1));
        end
        
        units_comb(end).ID = 0;
        units_comb(end).spikes = spikes;
    end
end

%

cds_new = [];

for f = fieldnames(cds)'
    cds_new.(f{1}) = cds.(f{1});
    
end
cds_new.units = units_comb;
cds = cds_new;