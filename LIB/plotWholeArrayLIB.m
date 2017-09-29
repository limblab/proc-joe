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
    
    plottedHere = zeros(10,10);
    pos = [0.1360,0.9040,0.0560,0.0560];
    
    for nn = 1:size(cds.units,2)
        % get grid spot according to optsArray.MAP_FILE
        mapIdx = find(cds.units(nn).chan == optsArray.ARRAY_MAP.chan);
        eRow = 11 - optsArray.ARRAY_MAP.row(mapIdx);
        subplotIdx = (eRow-1)*10 + optsArray.ARRAY_MAP.col(mapIdx);
        subplot(10,10,subplotIdx,'align')
        
        if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
            try
                feval(functionName,cds,nn,optsTask,optsPlot,optsSave);
                % remove axes
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                set(gca,'XLabel',[]);
                set(gca,'YLabel',[]);
                set(gca,'XTick',[]);
                set(gca,'XMinorTick','off')
                set(gca,'YTick',[]);
                set(gca,'YMinorTick','off')
                plottedHere(eRow,optsArray.ARRAY_MAP.col(mapIdx)) = 1;
            catch
                hold on
                ax = gca;
                XLIM = ax.XLim;
                YLIM = ax.YLim;
                stimDotX = (XLIM(2) + XLIM(1))/2;
                stimDotY = (YLIM(2) + YLIM(1))/2;
                plot(stimDotX,stimDotY,'k.','markersize',75)
                xlim(XLIM)
                ylim(YLIM)
                set(ax,'visible','off')
                set(ax,'color','none')
                plottedHere(eRow,optsArray.ARRAY_MAP.col(mapIdx)) = 1;
            end
        else
            
        end
    end
    
    for i = 1:10
        for j = 1:10
            if(~plottedHere(i,j))
                % get grid spot according to optsArray.MAP_FILE
                eRow = i;
                subplotIdx = (eRow-1)*10 + j;
                subplot(10,10,subplotIdx,'align')
                
                COLOR = 'k';
                
                ax = gca;
                XLIM = ax.XLim;
                YLIM = ax.YLim;
                stimDotX = (XLIM(2) + XLIM(1))/2;
                stimDotY = (YLIM(2) + YLIM(1))/2;
                plot(stimDotX,stimDotY,'.','color',COLOR,'markersize',75)
                xlim(XLIM)
                ylim(YLIM)
                set(ax,'visible','off')
                set(ax,'color','none')
            end
        end
    end
    
    % stim electrode
    posIdx=find(optsArray.ARRAY_MAP.chan==optsArray.STIM_ELEC);
    eRow=optsArray.ARRAY_MAP.row(posIdx);
    eRow = 11 - eRow;
    eCol=optsArray.ARRAY_MAP.col(posIdx);
    subplot(10,10,10*(eRow-1)+eCol,'align');
    hold on
    ax = gca;
    XLIM = ax.XLim;
    YLIM = ax.YLim;
    if(plottedHere(eRow,eCol)==1)
        stimBoxX=[XLIM(1),XLIM(2),XLIM(2),XLIM(1),XLIM(1)];
        stimBoxY=[YLIM(1),YLIM(1),YLIM(2),YLIM(2),YLIM(1)];
        plot(stimBoxX,stimBoxY,'mp-','linewidth',3)
    else
        stimDotX = (XLIM(2) + XLIM(1))/2;
        stimDotY = (YLIM(2) + YLIM(1))/2;
        plot(stimDotX,stimDotY,'m.','markersize',75)
        xlim(XLIM)
        ylim(YLIM)
        set(ax,'visible','off')
        set(ax,'color','none')
    end
    
    
end


function [optsArray] = configureOptionsArray(optsArrayInput)

    optsArray.MAP_FILE = '';
    optsArray.STIM_ELEC = '';
    
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