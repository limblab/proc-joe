function [figureHandle] = plotWholeArrayLIB(cds,functionName,optsTask,optsPlot,optsSave,optsArray)
%% wrapper class that generates a 10x10 grid of figures and plots on a group of them
% according to functionName. Also takes in some options for the plots
    
    optsArray = configureOptionsArray(optsArray);
    optsArray.ARRAY_MAP = loadMapFile(optsArray.MAP_FILE);
    
    figureHandle = figure();
    % no labels or titles
    optsPlot.X_LABEL = '';
    optsPlot.Y_LABEL = '';
    optsPlot.PLOT_TITLE = 0;
    optsPlot.LEGEND_STRING = '';
    for nn = 1:size(cds.units,2)
        if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
            % get grid spot according to optsArray.MAP_FILE
            mapIdx = find(cds.units(nn).chan == optsArray.ARRAY_MAP.chan);
            subplotIdx = (optsArray.ARRAY_MAP.row(mapIdx)-1)*10 + optsArray.ARRAY_MAP.col(mapIdx);
            subplot(10,10,subplotIdx)
            feval(functionName,cds,nn,optsTask,optsPlot,optsSave);
        end
    end

end


function [optsArray] = configureOptionsArray(optsArrayInput)

    optsArray.MAP_FILE = '';

    %% check if in optsArray and optsArrayInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsArrayInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(optsArray,inputFieldnames{fn}))
               optsArray.(inputFieldnames{fn}) = optsArrayInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end