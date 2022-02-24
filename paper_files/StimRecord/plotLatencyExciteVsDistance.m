function [outputData,figureHandles,FITS,gof] = plotLatencyExciteVsDistance(arrayData,mapFileName,opts)

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
    
    distancesAll = cell(size(arrayData{1}.spikeTrialTimes,2),1);
    latencyAll = cell(size(arrayData{1}.spikeTrialTimes,2),1);
    suspectUnits = [];
    
    for chan = 1:size(arrayData{1}.spikeTrialTimes,1)
        for wave = 1:size(arrayData{1}.spikeTrialTimes,2)
            %% get data 
            latency = zeros(size(arrayData,1),1);
            distances = zeros(size(arrayData,1),1);
            STIMCHAN_POS = [11-MAP_DATA.row(find(MAP_DATA.chan == arrayData{1,1}.CHAN_LIST(chan))), ...
                MAP_DATA.col(find(MAP_DATA.chan == arrayData{1,1}.CHAN_LIST(chan)))];
            
            for unit = 1:numel(arrayData)
                if(arrayData{unit}.isExcitatory{chan,wave})
                    
                    latency(unit,1) = arrayData{unit}.excitatoryLatency{chan,wave}(2); % 2 is the peak
%                     if(latency(unit) > 8)
%                         disp(num2str(unit))
%                     end
                    distances(unit,1) = 400*sqrt((arrayData{unit}.ROW-STIMCHAN_POS(1)).^2 + (arrayData{unit}.COL-STIMCHAN_POS(2)).^2);
                    
                    if(latency(unit,1) < 2 && distances(unit,1) > 2000)
                        suspectUnits(end+1,:) = [unit,chan,wave];
                    end
                end
            end
            latency(distances == 0) = [];
            distances(distances == 0) = [];
            
            distances(latency <= 0) = [];
            latency(latency <= 0) = [];
            
            distancesAll{wave} = [distancesAll{wave};distances];
            latencyAll{wave} = [latencyAll{wave};latency];
            
            %% plot distance vs latency
            if(~opts.PLOT_ON_ONE_FIGURE)
                figureHandles{end+1} = figure();
                c = 'k';
            else
                c = getColorFromList(1,wave-1);
            end
            
            plot(distances,latency,'.','markersize',opts.MARKER_SIZE,'color',c);

            if(~opts.PLOT_ON_ONE_FIGURE)
                f=gcf;
%                 xlim([0,f.Children.XLim(2)]);
%                 ylim([0,f.Children.YLim(2)])
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
    
    FITS = cell(numel(distancesAll),1); gof = cell(numel(distancesAll),1);
    for wave = 1:numel(distancesAll)
        [temp_fit,temp_gof] = fit(distancesAll{wave},latencyAll{wave},'a*x+b');
        FITS{wave} = temp_fit;
        gof{wave} = temp_gof;
            if(opts.PLOT_FIT)
                plot([0,3000],feval(FITS{wave},[0,3000]),'--','color',getColorFromList(1,wave),'linewidth',2);
            end
    end
    
    xlabel('Distance (\mum)');
    ylabel('Peak latency (ms)');
    set(gca,'fontsize',14)
    formatForLee(gcf)
    
    outputData.distancesAll = distancesAll;
    outputData.latencyAll = latencyAll;
    outputData.suspiciousUnits = unique(suspectUnits,'rows');
end


function [opts] = configureOpts(optsInput)

    
    opts.MARKER_SIZE = 20;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.PLOT_FIT = 1;
    
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