function [] = spindleStimEMGAnalysis(cds, eventTimes,filename)
% plot EMG for all muscles for all events (a single plot for each
% stimulation event)
% this is stored in cds.emg.t (for time stamps) and cds.emg.(muscle name)
% or cds.emg{:,colNum} to get data 

sampRate = 1/(cds.emg.t(2) - cds.emg.t(1));
for eventTime = eventTimes'
    % grab EMG data around eventTime -- say 500ms before and 500ms after
    timeBefore = 0.5; % in seconds
    timeAfter = 0.5; % in seconds
    [~,idxEvent] = min(abs(cds.emg.t-eventTime));
    idxEMG = [idxEvent - floor(timeBefore*sampRate), idxEvent+floor(timeAfter*sampRate)];
    emgData = cds.emg(idxEMG(1):idxEMG(2),:); % still a table
    
    % plot raw data?
%     [b,a] = butter(6,[53,67]/(sampRate/2),'stop'); % apparently it needs to be filtered
    figure;
    for i = 2:numel(emgData.Properties.VariableNames)
        subplot(numel(emgData.Properties.VariableNames)-1,1,i-1)
%         emgDataPlot = filtfilt(b,a,emgData{:,i});
        emgDataPlot = emgData{:,i};
        plot(emgData{:,1},emgDataPlot)
    end
end

end