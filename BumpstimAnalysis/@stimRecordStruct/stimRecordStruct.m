classdef stimRecordStruct < matlab.mixin.SetGet & operationLogger
    
    
    properties (SetAccess = private, GetAccess = public)
        meta
        units
        trials
        kin
        stimData
    end
    
    methods (Static = true)
        function srs = stimRecordStruct(varargin) % constructor
            % meta related data
            set(srs,'meta',[]);
                        
            % empty units
            srs.units = {};
            srs.kin = {};
            srs.trials = {};
            srs.stimData = {};
        end
    end
    
    methods
        % set methods -- really not meant for use once everything is moved
        % over
    
        function set.meta(srs,meta)
            srs.meta = meta;
        end
        function set.kin(srs,kin)
            srs.kin = kin;
        end
        function set.trials(srs,trials)
            srs.trials = trials;
        end
        function set.stimData(srs,stimData)
            srs.stimData = stimData;
        end
    end
    
    methods
        % append to data -- this is what should be used
        function cdsToStimRecordStruct(srs,cds,opts)
            % this function converts a given cds to a stimRecordStruct
            % this copies trials, kin and meta information over
            % this also adds stim meta information to meta.stim
            % this also removes unsorted or invalid units
            
            % copy kin and trial info over
            set(srs,'kin',cds.kin);
            set(srs,'trials',cds.trials);
            
            % copy meta and set meta.stim
            metaTemp = cds.meta;
            metaTemp.stim.waveformInformation = cds.waveforms;
            metaTemp.stimOn = cds.stimOn;
            set(srs,'meta',metaTemp);
            
            % copy artifact data over
        end
        
    end
    
    methods (Static = false)
        % methods are located elsewhere -- must be in this folder
        addData(srs,folderpath,filenames,varargin)
    end
end

