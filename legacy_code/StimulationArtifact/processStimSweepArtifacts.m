%% process stim sweep artifacts

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\ArtificialMonkey_20171221\Chips\';
inputData.mapFile = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

inputData.monkey = 'monkeyChips';
inputData.lab = 6;
inputData.ranBy = 'ranByJoe';
inputData.array1 = 'arrayLeftS1';
inputData.task = 'taskWF';

inputData.badChList = [];
inputData.interpulse = 0.000053;
inputData.pWidth1 = 0.0002;
inputData.pWidth2 = 0.0002;

inputData.windowSize=30*99;%in points
inputData.presample=60;%in points
inputData.plotRange=0.8;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];


pwd=cd;
% fileList = dir('*Chips*chan10*.nev');

%% run process script to get artifact data
[~,outputData] = processStimArtifact(folderpath,inputData);
cd(pwd);
%% determine length of 'hitting the rails' after 2ms on average for each amplitude
railValue = 8000; % a little below the true rail 
railOffset = 2*30; % ignore the first 2ms

for art = 1:numel(outputData.artifactData)
    for amp = 1:100
        artifactIdx = amp:100:size(outputData.artifactData(art).artifact,2);
        stimChannelIdx = find(outputData.chList == outputData.artifactData(art).stimChannel);
        artifact = squeeze(outputData.artifactData(art).artifact(stimChannelIdx,artifactIdx,:));
        
        for a = 1:size(artifact,1)  
            railIdx = find(abs(artifact(a,inputData.presample+railOffset:end)) < railValue);
            outputData.artifactData(art).railIdxs(a) = railIdx(1);
        end
        outputData.artifactData(art).timeOnRails(amp,1) = (floor(mean(outputData.artifactData(art).railIdxs)) + railOffset)/30;
    end
end


