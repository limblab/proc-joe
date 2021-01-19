% combine outputData's so that there are more examples
% load in inputData and provide folderpath

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Duncan\MultiElec\Duncan_20181211_multiElec\';
fpath = cd;
cd(folderpath)
file_list = dir('*outputData*');
cd(fpath)
outputData_all = [];

for f = 1:10%numel(file_list)
    disp(f)
    load([file_list(f).folder, '\', file_list(f).name]);
    
    if(f == 1)
        outputData_all = outputData;
    else
        % artifact data
        outputData_all.artifactData.t = [outputData_all.artifactData.t; outputData.artifactData.t];
        outputData_all.artifactData.artifact = [outputData_all.artifactData.artifact; outputData.artifactData.artifact];
        % waveforms
        outputData_all.waveforms.waveSent = [outputData_all.waveforms.waveSent; outputData.waveforms.waveSent];
        outputData_all.waveforms.chanSent = [outputData_all.waveforms.chanSent; outputData.waveforms.chanSent];
    end
end
