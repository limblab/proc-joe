function [heatmapPD,PDscaled] = plotHeatmapsPD(td_all,pd_all,mapData,optsPD)

    % configure opts and set default values
    optsPD = configureOptsPD(optsPD);
    
    % get heatmap data
    heatmapDataPD = {};
    [heatmapDataPD,alphaArray,stimRow,stimCol,PDscaled] = getHeatmapDataPD(td_all,pd_all,optsPD,mapData);
    
    %plot the pds
    f = figure();
    f.Name = strcat(optsPD.FIGURE_PREFIX,'_stimChan',num2str(optsPD.CENTER_CHANNEL),'_preferredDirectionsHeatmap');
    heatmapPD = imagesc(heatmapDataPD,'alphaData',alphaArray);
    
    colormap(flip(viridis,1));
    colorbar;
    
    % magenta box for stim chan
    hold on
    rectangle('position',[(stimCol-0.5) (stimRow-0.5) 1 1],'edgecolor','m','linewidth',2);
    plot([stimCol-0.5,stimCol+0.5],[stimRow-0.5,stimRow+0.5],'m','linewidth',2);
    plot([stimCol+0.5,stimCol-0.5],[stimRow-0.5,stimRow+0.5],'m','linewidth',2);
    
end

function [heatmapDataPD,alphaArray,stimRow,stimCol,PDscaled] = getHeatmapDataPD(td_all,pd_all,optsPD,mapData)

    %create array that you're going to pass to imagesc
    heatmapDataPD = nan(optsPD.NUM_ROWS,optsPD.NUM_COLS);
    alphaArray = zeros(optsPD.NUM_ROWS,optsPD.NUM_COLS);

    %ensure all angles are between 1-180deg
    PDtemp = radtodeg(pd_all.velPD);
    
    %calculate pd's relative to stimulated channel pd (need to change this
    %to hardcoded input angle)
    for i=1:numel(PDtemp)
        PDtemp(i) = angleDiff(PDtemp(optsPD.CENTER_CHANNEL),PDtemp(i),false,false);
    end 

    %put the angles in their respective locations in the array
    for pd=optsPD.PLOT_CHANNELS
        for chan=1:96
            if td_all(1).LeftS1_unit_guide(pd) == mapData.chan(chan)
                heatmapDataPD((11-mapData.row(chan)),mapData.col(chan)) = PDtemp(pd);
                alphaArray(11-mapData.row(chan),mapData.col(chan))=1;
                %setting stimrow and stimcol
                if mapData.chan(chan) == optsPD.CENTER_CHANNEL
                    stimRow = 11 - mapData.row(chan);
                    stimCol = mapData.col(chan);
                end
            end
        end
    end
    
    %scale from -1 to 1
    heatmapDataPDscaled = heatmapDataPD;
    heatmapDataPDscaled = -2.*(heatmapDataPDscaled-(min(heatmapDataPDscaled)))./(max(heatmapDataPDscaled)-min(heatmapDataPDscaled)) + 1;
    
    %save variable for combineHeatmaps
    PDscaled = heatmapDataPDscaled;

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