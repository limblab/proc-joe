<<<<<<< HEAD
function [ filteredlfp ] = filterCDSData( cds )
% applies a filter to the cds data in lfp, then returns the cds
% this is a preliminary function for the stim artifact project to see if
% what I am doing is right

sampRate = 1./(cds.lfp{2:100,1}-cds.lfp{1:99,1});
sampRate = mean(sampRate);
fc = 10;
[b,a] = butter(2,fc/(sampRate/2),'low');
% apply filter to all electrodes
filteredlfp = cds.lfp;
for i = 2:size(cds.lfp,2)
    % reverse the signal
    signalToFilter = fliplr(cds.lfp{:,i}')';
    % filter and unreverse signal
%     filteredlfp{:,i} = fliplr(filter(b,a,signalToFilter)')';
    filteredlfp{:,i} = fliplr(signalToFilter')';
end


end

=======
function [ filteredlfp ] = filterCDSData( cds )
% applies a filter to the cds data in lfp, then returns the cds
% this is a preliminary function for the stim artifact project to see if
% what I am doing is right

sampRate = 1./(cds.lfp{2:100,1}-cds.lfp{1:99,1});
sampRate = mean(sampRate);
fc = 10;
[b,a] = butter(2,fc/(sampRate/2),'low');
% apply filter to all electrodes
filteredlfp = cds.lfp;
for i = 2:size(cds.lfp,2)
    % reverse the signal
    signalToFilter = fliplr(cds.lfp{:,i}')';
    % filter and unreverse signal
%     filteredlfp{:,i} = fliplr(filter(b,a,signalToFilter)')';
    filteredlfp{:,i} = fliplr(signalToFilter')';
end


end

>>>>>>> 0e441d3b4e9a2e8831546d9f89425fae8b961b8d
