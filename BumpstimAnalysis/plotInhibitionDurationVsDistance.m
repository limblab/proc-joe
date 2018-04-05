function [figureHandles,FITS] = plotInhibitionDurationVsDistance(arrayData,mapFileName,opts)
    %% configure opts and set default values
    opts = configureOpts(opts);
    
    %% useful constants
    figureHandles = {};
    FITS = [];
    MAP_DATA = loadMapFile(mapFileName);

    if(opts.PLOT_ON_ONE_FIGURE)
        figureHandles{end+1} = figure();
        hold on
    end
    
    distancesAll = [];
    inhibitionDurationAll = [];
    
    for chan = 1:size(arrayData{1}.spikeTrialTimes,1)
        for wave = 1:size(arrayData{1}.spikeTrialTimes,2)
            %% get data 
            inhibitionDuration = zeros(size(arrayData,1),1);
            distances = zeros(size(arrayData,1),1);
            STIMCHAN_POS = [MAP_DATA.row(find(MAP_DATA.chan == arrayData{1,1}.CHAN_LIST(chan))), MAP_DATA.col(find(MAP_DATA.chan == arrayData{1,1}.CHAN_LIST(chan)))];
            
            for unit = 1:size(arrayData,1)
                inhibitionDuration(unit) = diff(arrayData{unit}.inhibitoryLatency{chan,wave});
                distances(unit) = 400*sqrt((arrayData{unit,1}.ROW-STIMCHAN_POS(1)).^2 + (arrayData{unit,1}.COL-STIMCHAN_POS(2)).^2);
            end
            distances(inhibitionDuration <= 0) = [];
            inhibitionDuration(inhibitionDuration <= 0) = [];
            
            distancesAll = [distancesAll;distances];
            inhibitionDurationAll = [inhibitionDurationAll;inhibitionDuration];
            
            %% plot distance vs inhibition duration
            if(~opts.PLOT_ON_ONE_FIGURE)
                figureHandles{end+1} = figure();
                c = 'k';
            else
                c = opts.COLORS{chan*wave};
            end
            plot(distances,inhibitionDuration,'.','markersize',opts.MARKER_SIZE,'color',c);

            
            %% save figures
            if(opts.FIGURE_SAVE && strcmpi(FIGURE_DIR,'')~=1)
                FIGURE_NAME = strcat(FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),'_heatmap');
                saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,FIGURE_NAME);
            end
            
            
        end
    end
    
    FITS = fit(distancesAll,inhibitionDurationAll,'a*x+b');
end


function [opts] = configureOpts(optsInput)

    
    opts.MARKER_SIZE = 20;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.PLOT_ON_ONE_FIGURE = 1;
    opts.COLORS = {'r','b','g','k','m',[0.5,0.5,0.5],'y'};

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