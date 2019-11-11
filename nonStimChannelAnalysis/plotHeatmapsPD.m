function [heatmapPD,heatmapDataPD,alphaArray] = plotHeatmapsPD(td_all,pd_all,mapData,optsPD)

    % configure opts and set default values
    optsPD = configureOptsPD(optsPD);
    
    % get heatmap data in 10x10 array
    heatmapDataPD = {};
    [heatmapDataPD,~] = getHeatmapDataPD(td_all,pd_all,mapData);
    
    % only plot entries in optsPD.PLOT_CHANNELS by setting alphaArray() =
    % 1;
    alphaArray = zeros(size(heatmapDataPD));
    
    for plotChan = 1:numel(optsPD.PLOT_CHANNELS)
        mapDataIdx = find(mapData.chan == optsPD.PLOT_CHANNELS(plotChan));
        alphaArray(11-mapData.row(mapDataIdx),mapData.col(mapDataIdx)) = 1;
    end
    
    %plot the pds
    f = figure();
    f.Name = strcat(optsPD.FIGURE_PREFIX,'_preferredDirectionsHeatmap');
    heatmapPD = imagesc(heatmapDataPD,'alphaData',alphaArray);
    axis square
    colormap(colorcet('C9'));
    colorbar;
    
%     % magenta box for stim chan
%     hold on
%     rectangle('position',[(stimCol-0.5) (stimRow-0.5) 1 1],'edgecolor','m','linewidth',2);
%     plot([stimCol-0.5,stimCol+0.5],[stimRow-0.5,stimRow+0.5],'m','linewidth',2);
%     plot([stimCol+0.5,stimCol-0.5],[stimRow-0.5,stimRow+0.5],'m','linewidth',2);
    
end


function [optsPD] = configureOptsPD(optsPDInput)
    optsPD.MAKE_BAR_PLOT = 1;

    optsPD.PLOT_CHANNELS = [1:96];
    optsPD.CENTER_CHANNEL = 1;

    optsPD.MAX_RATIO = 1;
    optsPD.MIN_RATIO = -0.2;
    optsPD.LOG_SCALE = 0;
    optsPD.LOG_PARAM = 9;
    
    optsPD.FIGURE_SAVE = 0;
    optsPD.FIGURE_DIR = '';
    optsPD.FIGURE_PREFIX = '';
    
    optsPD.NUM_ROWS = 10;
    optsPD.NUM_COLS = 10;
    
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsPDInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(optsPD,inputFieldnames{fn}))
               optsPD.(inputFieldnames{fn}) = optsPDInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end