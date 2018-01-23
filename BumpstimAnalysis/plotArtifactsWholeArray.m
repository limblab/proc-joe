function [figHandles] = plotArtifactsWholeArray(cds,mapFilename, opts)
% plot artifacts across the array (or a subset)

    opts = configureOptions(opts);

    % make figure
    figure('Position',[357.8000 56.2000 838.4000 726.4000]);
    mapData=loadMapFile(mapFilename);
    %establish trackers for information that is common to all files:
    posList = [mapData.row,mapData.col];
    plottedHere = zeros(10,10);% plotted neurons here
    
    numRows = max(opts.ROWS_PLOT) - min(opts.ROWS_PLOT) + 1;
    numCols = max(opts.COLS_PLOT) - min(opts.COLS_PLOT) + 1;
    
    % get artifact indexes corresponding to waveform to plot
    artifactMult = numel(cds.stimOn)/size(cds.artifactData.artifact,1);
    channelList = unique(cds.waveforms.chanSent);
    artifactIdx = find(cds.waveforms.waveSent(1:artifactMult:end) == opts.WAVEFORM_PLOT & ...
        cds.waveforms.chanSent(1:artifactMult:end) == channelList(opts.CHANNEL_PLOT));
    
    for chan = 1:96
        posIdx=find(mapData.chan==chan);
        eRow=posList(posIdx,1);
        eRow = 11 - eRow;
        eCol=posList(posIdx,2);
        if(~isempty(find(opts.ROWS_PLOT,eRow)) && ~isempty(find(opts.COLS_PLOT,eCol)))
            h=subplot(numRows,numCols,numCols*(eRow-min(opts.ROWS_PLOT))+eCol);
            plottedHere(eRow,eCol) = 1;
            
            xData = ((0:1:size(cds.artifactData.artifact,3)-1))/30;
            
            if(opts.MARK_BANK)
                switch mapData.bank{posIdx}
                    case 'A'
                        colorPlot = 'r';
                    case 'B' 
                        colorPlot = [0, 0.5, 0];
                    case 'C'
                        colorPlot = 'b';
                end
            else
                colorPlot = 'r';
            end
            plot(xData,squeeze(cds.artifactData.artifact(artifactIdx(1:min(opts.ARTIFACTS_PLOT,numel(artifactIdx))),chan,:))','color',colorPlot,'linewidth',1.5)
        
            % remove axes
            xlim(opts.XLIM)
            ylim(opts.YLIM)
            
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            set(gca,'XLabel',[]);
            set(gca,'YLabel',[]);
            set(gca,'XTick',[]);
            set(gca,'XMinorTick','off')
            set(gca,'YTick',[]);
            set(gca,'YMinorTick','off')
            set(gca,'Visible','off')
        end

    end

    % mark stim chan (if requested)
    if(opts.MARK_STIM_CHAN)
        posIdx=find(mapData.chan==channelList(opts.CHANNEL_PLOT));
        eRow=posList(posIdx,1);
        eRow = 11 - eRow;
        eCol=posList(posIdx,2);
        h=subplot(10,10,10*(eRow-1)+eCol);
        % h=subplot(5,3,3*(eRow-4)+eCol);
        hold on
        ax = h;
        XLIM = ax.XLim;
        YLIM = ax.YLim;
        if(plottedHere(eRow,eCol)==1)
            stimBoxX=[XLIM(1),XLIM(2),XLIM(2),XLIM(1),XLIM(1)];
            stimBoxY=[YLIM(1),YLIM(1),YLIM(2),YLIM(2),YLIM(1)];
            plot(ax,stimBoxX,stimBoxY,'m-','linewidth',2)
        else
        %                 stimCrossX1 = [XLIM(1),XLIM(2)];
        %                 stimCrossY1 = [YLIM(1),YLIM(2)];
        %                 stimCrossX2 = [XLIM(1),XLIM(2)];
        %                 stimCrossY2 = [YLIM(2),YLIM(1)];
            stimDotX = (XLIM(2) + XLIM(1))/2;
            stimDotY = (YLIM(2) + YLIM(1))/2;
            plot(ax,stimDotX,stimDotY,'m','markersize',75)
            xlim(XLIM)
            ylim(YLIM)
            set(ax,'visible','off')
            set(ax,'color','none')
        end
    end

end

function [opts] = configureOptions(optsInput)
    opts.ROWS_PLOT = 1:10;
    opts.COLS_PLOT = 1:10;
    opts.CHANNEL_PLOT = 1;
    opts.WAVEFORM_PLOT = 1;
    opts.ARTIFACTS_PLOT = 10;
    
    opts.XLIM = [-0.3,5];
    opts.YLIM = [-9000, 9000];
    opts.MARK_STIM_CHAN = 1;
    opts.MARK_BANK = 1;
    %% check if in opts and optsInput, overwrite if so
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


