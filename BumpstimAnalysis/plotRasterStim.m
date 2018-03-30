function [figureHandle] = plotRasterStim(unitData,NEURON_NUMBER,optsPlot)

    %% configure opts and set default values
    opts = configureOpts(optsPlot);
    
    NUM_WAVEFORM_TYPES = size(unitData.spikeTrialTimes,2);
    NUM_CHANS = size(unitData.spikeTrialTimes,1);
    CHAN_LIST = unitData.CHAN_LIST;
    

    %% format data and call general plot function
    
    for chan = 1:NUM_CHANS
        for wave = 1:NUM_WAVEFORM_TYPES
            
            xData = unitData.spikeTrialTimes{chan,wave}*1000;
            yData = unitData.stimData{chan,wave};

            if(strcmpi(opts.SORT_DATA,'velocity')) % handle the sort
                vx = unitData.kin{chan,wave}.vx;
                vy = unitData.kin{chan,wave}.vy;
                vel = mean(vx.^2 + vy.^2,2);
                [~,sortIdx] = sort(vel);
                [xData,yData] = sortRasterData(xData,yData,sortIdx);
            elseif(strcmpi(opts.SORT_DATA,'acceleration'))
                ax = unitData.kin{chan,wave}.ax;
                ay = unitData.kin{chan,wave}.ay;
                accel = mean(ax.^2 + ay.^2,2);
                [~,sortIdx] = sort(accel);
                [xData,yData] = sortRasterData(xData,yData,sortIdx);
            else % plotRasterLIB handles the sort otherwise
                optsPlot.SORT_DATA = opts.SORT_DATA;
            end
            % plot related information
            optsPlot.MAKE_FIGURE = opts.MAKE_FIGURE;
            optsPlot.X_LIMITS = [-opts.PRE_TIME*1000,opts.POST_TIME*1000];
            optsPlot.Y_LIMITS = [-3,unitData.numStims(chan,wave)+1];
            if(~opts.MAKE_SUBPLOTS || wave == NUM_WAVEFORM_TYPES)
                optsPlot.X_LABEL = 'Time after stimulation onset (ms)';
            else
                optsPlot.X_LABEL = '';
            end
            optsPlot.Y_LABEL = 'Stimuli';
            % deals with title requests
            if(opts.PLOT_TITLE)
                if(strcmp(opts.TITLE_TO_PLOT,'') == 0)
                    optsPlot.TITLE = titleToPlot;
                else
                    optsPlot.TITLE = strcat('Stim Chan: ',num2str(CHAN_LIST(chan)),' Wave: ',num2str(wave));
                end
            end
            optsPlot.Y_TICK = [0;max(unitData.numStims(chan,wave))];
            optsPlot.Y_MINOR_TICK = 'off';
            optsPlot.Y_TICK_LABEL = {num2str(1),num2str(unitData.numStims(chan,wave))};
            if(strcmpi(opts.SORT_DATA,'velocity'))
                optsPlot.Y_TICK_LABEL = {'Low','High'};
                optsPlot.Y_LABEL = 'Velocity';
            end
            optsSave.FIGURE_SAVE = opts.FIGURE_SAVE;
            optsSave.FIGURE_DIR = opts.FIGURE_DIR;
            if(opts.FIGURE_SAVE)
%                 optsSave.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'nn',num2str(NEURON_NUMBER),'_chan',num2str(unitData.CHAN_REC),'_stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'_raster');
                optsSave.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'nn',num2str(NEURON_NUMBER),'_chan',num2str(unitData.CHAN_REC),'_waveNum',num2str(wave),'_raster');
            else
                optsSave.FIGURE_NAME = '';
            end

            optsPlot.MARKER_STYLE = opts.MARKER_STYLE;
            optsPlot.MARKER_SIZE = opts.MARKER_SIZE;
            
            % plot Raster using general function
            figureHandle{chan,wave} = plotRasterLIB(xData,yData,optsPlot,optsSave);
   
        end
    end
    



end

function [xDataSorted,yDataSorted] = sortRasterData(xData,yData,sortIdx)

    yDataSorted = yData;
    xDataSorted = xData;
    counter = 1;
    for sIdx = 1:numel(sortIdx)
        yDataSorted(counter:counter+sum(yData==sortIdx(sIdx))-1) = sIdx;
        xDataSorted(counter:counter+sum(yData==sortIdx(sIdx))-1) = xData(yData==sortIdx(sIdx));
        counter = counter+sum(yData==sortIdx(sIdx));
    end
    
end


function [opts] = configureOpts(optsInput)

    opts = [];
    opts.MARKER_SIZE = 4;
    opts.MARKER_STYLE  = '.';
    opts.PLOT_TITLE = 1;
    opts.TITLE_TO_PLOT = '';
    opts.MAKE_FIGURE = 1;
    opts.SORT_DATA = 0;
    opts.MAKE_SUBPLOTS = 0;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.PRE_TIME = 5/1000;
    opts.POST_TIME = 30/1000;
    
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end


% 
% %% plot raw waveforms around the artifact
%     if(plotSpikeWaveforms == 1)
%         for chan = chansPlot
%             for fig = waveformTypesPlot
%                 plotWaveformsStim(cds,NEURON_NUMBER,chan,fig,'timeAfterStimRawNoStim',timeAfterStimRawNoStim,'timeAfterStimRawArtifact',timeAfterStimRawArtifact,...
%                     'makeFigure',1,'plotTitle',plotTitle,'title',titleToPlot,'stimElectrode',stimElectrode,'saveFigures',saveFigures,...
%                     'figDir',figDir,'figPrefix',figPrefix,'plotFiltered',plotWaveformFiltered,'maxWaveformsPlot',maxWavesPlot);
%             end
%         end
%     end
% 
%     %% plot sample of artifacts with waveform and without waveform raw
%     if(plotArtifacts)
%         for pAF = plotArtifactFiltered
%             if(pAF == 0 || pAF == 1)
%                 for chan = chansPlot
%                     for fig = waveformTypesPlot
%                         plotArtifactsStim(cds,NEURON_NUMBER,chan,fig,'plotTitle',plotTitle,'title',titleToPlot,...
%                             'maxArtifactsPerPlot',maxArtifactsPerPlot,'timeAfterStim',timeAfterStimRawArtifact,...
%                             'rowSubplot',rowSubplotArtifact,'colSubplot',colSubplotArtifact,'plotFiltered',pAF,...
%                             'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,'templateSubtract',0,...
%                             'randomSample',0);
%                     end
%                 end
%             end
%         end
%     end
