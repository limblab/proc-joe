function [outputDataFiltered] = performFilteringStep( outputData, inputData, highPassCutoff, lowPassCutoff, filterOrder, folderpath )

sampRate = 30000; % hz
fhigh = highPassCutoff; % hz
flow = lowPassCutoff; % hz
disp(strcat('high',num2str(fhigh),'low',num2str(flow),'order',num2str(filterOrder)))
doNothing = 0;
noFilter = 0;

% build filter according to cutoffs provided
if(flow==-1 && fhigh == -1) % no filter case
    doNothing = 0;
    noFilter = 1;
elseif(fhigh == -1) % only use low pass
    [filterStruct.b,filterStruct.a] = butter(filterOrder,flow/(sampRate/2),'low');
elseif(flow == -1) % only use high pass
    [filterStruct.b,filterStruct.a] = butter(filterOrder,fhigh/(sampRate/2),'high');
elseif(flow <= fhigh) % do nothing because invalid filter
    doNothing = 1;
else % band pass
    [filterStruct.b,filterStruct.a] = butter(filterOrder,[fhigh,flow]/(sampRate/2),'bandpass');
end

if(~doNothing && ~noFilter) % filter and save png and .mat file
    outputDataFiltered = filterArtifactData(outputData,'filter',filterStruct);
    % plot results and save png
%     plotArtifactsAllStimChannels(outputDataFiltered,inputData,folderpath,'Name',filterStruct.name,'noPlots',1);
%     save(strcat(folderpath,filesep,'Raw_Figures',filesep,filterStruct.name(2:end)),'outputDataFiltered','inputData');
    close all
elseif(~doNothing && noFilter) % no filter, save .mat file
%     plotArtifactsAllStimChannels(outputData,inputData,folderpath,'Name','_noFilter','noPlots',1);
    outputDataFiltered = outputData;
%     save(strcat(folderpath,filesep,'Raw_Figures',filesep,'noFilter'),'outputDataFiltered','inputData');
    close all
end          

% % make sure to remove bad channels from summary statistic
% meanSettlingTime = mean(settlingTimeMat,3);
% meanSettlingTime = mean(meanSettlingTime(:,33:end),2);
% [meanSettlingTime, idx] = sort(meanSettlingTime);
% filterParamsMat = filterParamsMat(idx,:);
% % put everything in a structure to return
% filterResults.filters = filterParamsMat;
% filterResults.meanSettlingTime = [meanSettlingTime,filterParamsMat];
% filterResults.settlingTimes = settlingTimeMat(idx,:,:);
end

