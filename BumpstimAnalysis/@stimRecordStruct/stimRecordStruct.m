classdef stimRecordStruct < matlab.mixin.SetGet & operationLogger
    
    
    properties (SetAccess = public, GetAccess = public)
        monkey
        array
        dates
        numStims
        stimParameters
        units
    end
    
    methods (Static = true)
        function srs = stimRecordStruct(varargin)
            % meta related data
            d = {};
            monk = '';
            arr = '';
            nS = [];
            set(srs,'date',d);
            set(srs,'monkey',monk);
            set(srs,'array',arr);
            set(srs,'numStims',nS);
            
            % empty stimParameters
            srs.stimParameters = cell2table(cell(0,7),'VariableNames',{'stimChan','stimElec','amp','pw1','pw2','interphase','interpulse'});
            
            % empty units
            srs.units = {};
        end
    end
    
    methods
        % set methods -- really not meant for use other than the meta field
        function set.stimParameters(srs,stimParameters)
            srs.stimParameters = stimParameters;
        end
        
        function set.units(srs,units)
            srs.units = units;
        end
        function set.monkey(srs,monkey)
            srs.monkey = monkey;
        end
        function set.array(srs,array)
            srs.array = array;
        end
        function set.dates(srs,dates)
            srs.dates = dates;
        end
        function set.numStims(srs,numStims)
            srs.numStims = numStims;
        end
        

    end
    
    methods
        % append to data -- this is what should be used
        function appendStimParameters(srs,stimParameters)
            if(isempty(srs.stimParameters))
                srs.stimParameters = cell2table(stimParameters,'VariableNames',{'chan','elec','amp','pw1','pw2','interphase','interpulse'});
            else
                srs.stimParameters = [srs.stimParameters; stimParameters];
            end
        end
        function appendUnits(srs,units)
            srs.units{end+1,1} = units;
        end
        function appendDate(srs,date)
            srs.dates{end+1,1} = date;
        end
        function appendNumStims(srs,numStims)
            srs.numStims(end+1,1) = numStims;
        end
    end
    
    methods (Static = false)
        % methods are located elsewhere -- must be in this folder
        addData(srs,folderpath,filenames,varargin)
    end
end

