function [] = addData( srs,folderpath,filenames,varargin )
% this function adds data to the stim record structure
% folderpath points to a folder, filenames point to the files to use --
% there can be more than 1 file in a folder

pwd = cd;
cd(folderpath)
if(~iscell(filenames))
    filenames = {filenames};
end
for fname = filenames
    fname = fname{1};
    if(exist(fname)~=0)
        load(fname)
        chanList  = unique(cds.waveforms.chanSent);
        for chan = 1:numel(unique(cds.waveforms.chanSent))
            for wave = 1:numel(unique(cds.waveforms.waveSent))
                % append information about file
                srs.appendDate(cds.meta.dateTime(1:strfind(cds.meta.dateTime,' ')-1))
                labelIdx = -1;
                for l = 1:size(cds.units,2)
                    if(cds.units(l).chan == chanList(chan))
                        labelIdx = l;
                    end
                end
                if(labelIdx~=-1)
                    stimParameters = {chanList(chan),cds.units(labelIdx).label,cds.waveforms.parameters(wave).amp1,cds.waveforms.parameters(wave).pWidth1,...
                        cds.waveforms.parameters(wave).pWidth2,cds.waveforms.parameters(wave).interphase,cds.waveforms.parameters(wave).interpulse};
                else
                    stimParameters = {chanList(chan),'',cds.waveforms.parameters(wave).amp1,cds.waveforms.parameters(wave).pWidth1,...
                        cds.waveforms.parameters(wave).pWidth2,cds.waveforms.parameters(wave).interphase,cds.waveforms.parameters(wave).interpulse};
                end
                srs.appendStimParameters(stimParameters);
                srs.appendNumStims(sum(cds.waveforms.waveSent == wave & cds.waveforms.chanSent == chanList(chan)));
                units = cds.units;
                
                % remove unsorted and invalid units
                removeList = [];
                for u = 1:size(units,2)
                    if(units(u).ID == 0 || units(u).ID == 255)
                        removeList(end+1,1) = u;
                    end
                end
                units(removeList) = [];
                % extract spikes that correspond to stimuli only -- also
                % orient based on stimuli
                timeAfterStimuli = min(diff(cds.stimOn)) - 20/1000;
                timeBeforeStimuli = 20/1000;
                for u = 1:size(units,2)
                    spikes = [];
                    stimNum = [];
                    for st = 1:numel(cds.stimOn)
                        if(cds.waveforms.chanSent(st) == chanList(chan) && cds.waveforms.waveSent(st) == wave)
                            mask = units(u).spikes.ts > cds.stimOn(st)-timeBeforeStimuli & units(u).spikes.ts < cds.stimOn(st)+timeAfterStimuli;
                            spikes = [spikes;units(u).spikes.ts(mask)-cds.stimOn(st)];
                            stimNum = [stimNum;st*ones(sum(mask),1)];
                        end
                    end
                    units(u).spikes = spikes;
                    units(u).stimNum = stimNum;
                end
                
                srs.appendUnits(units);
            end
        end
        
        
    else
        warning(strcat('File ', fname,' does not appear to exist. File will be skipped'));
    end
end

cd(pwd)
end

