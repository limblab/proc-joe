function [figureHandles,FITS] = plotLatencyExciteVsDistance(arrayData,mapFileName,opts)

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
    latencyAll = [];
    
    for chan = 1:size(arrayData{1}.spikeTrialTimes,1)
        for wave = 1:size(arrayData{1}.spikeTrialTimes,2)
            %% get data 
            latency = zeros(size(arrayData,1),1);
            distances = zeros(size(arrayData,1),1);
            STIMCHAN_POS = [11-MAP_DATA.row(find(MAP_DATA.chan == arrayData{1,1}.CHAN_LIST(chan))), MAP_DATA.col(find(MAP_DATA.chan == arrayData{1,1}.CHAN_LIST(chan)))];
            
            for unit = 1:size(arrayData,1)
                if(arrayData{unit}.isExcitatory{chan,wave})
                    latency(unit) = arrayData{unit}.excitatoryLatency{chan,wave}(2); % 2 is the peak
                    if(latency(unit) > 8)
                        disp(num2str(unit))
                    end
                end
                distances(unit) = 400*sqrt((arrayData{unit,1}.ROW-STIMCHAN_POS(1)).^2 + (arrayData{unit,1}.COL-STIMCHAN_POS(2)).^2);
            end
            distances(latency <= 0) = [];
            latency(latency <= 0) = [];
            
            distancesAll = [distancesAll;distances];
            latencyAll = [latencyAll;latency];
            
            %% plot distance vs latency
            if(~opts.PLOT_ON_ONE_FIGURE)
                figureHandles{end+1} = figure();
                c = 'k';
            else
                c = opts.COLORS{chan*wave};
            end
            plot(distances,latency,'.','markersize',opts.MARKER_SIZE,'color',c);

            if(~opts.PLOT_ON_ONE_FIGURE)
                f=gcf;
                xlim([0,f.Children.XLim(2)]);
                ylim([0,f.Children.YLim(2)])
                %% save figures
                if(opts.FIGURE_SAVE && strcmpi(FIGURE_DIR,'')~=1)
                    FIGURE_NAME = strcat(FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),'_LatencyExciteVsDistance');
                    saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,FIGURE_NAME);
                end
            end
            
            
        end
    end
    
    if(opts.PLOT_ON_ONE_FIGURE)
        f=gcf;
        xlim([0,f.Children.XLim(2)]);
        ylim([0,f.Children.YLim(2)]);
        %% save figures
        if(opts.FIGURE_SAVE && strcmpi(FIGURE_DIR,'')~=1)
            FIGURE_NAME = strcat(FIGURE_PREFIX,'_LatencyExciteVsDistance');
            saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,FIGURE_NAME);
        end
    end
    
    FITS = fit(distancesAll,latencyAll,'a*x+b');
end


function [opts] = configureOpts(optsInput)

    
    opts.MARKER_SIZE = 20;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.PLOT_ON_ONE_FIGURE = 1;
    opts.COLORS = {'k','r','b',[0,0.5,0],'m'};

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