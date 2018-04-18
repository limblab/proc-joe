function [ figureHandle ] = plotArrayMap( cds,NEURON_NUMBER,MAP_FILE_NAME, opts )
% this function plots things on a 10x10 grid of the array. Namely the
% stimulated and recording channel, though other functionality may be added
% in the future

    %% configure opts and set default values
    opts = configureOpts(opts);

    %% load map file
    arrayMap=loadMapFile(MAP_FILE_NAME);

    figureHandle = figure();
    %% remove axis labels
    ax = gca;
    ax.YTickLabel = {};
    ax.XTickLabel = {};

    %% create ROWSxCOLS grid
    surface(zeros(opts.NUM_ROWS+1,opts.NUM_COLS+1))
    colormap(opts.COLOR_MAP); 
    %% make square
    figureHandle.Position(4) = figureHandle.Position(3);
    figureHandle.Position(2) = figureHandle.Position(2) - 200; % move down to not be annoyingly off my screen
    
    ax = gca;
    ax.Children(1).LineWidth = 1.5; % thicker box boundaries
    %% label stim electrode -- currently using an 'S' though this can be changed maybe later
    if(~isempty(opts.STIM_ELECTRODE))
        for stElec = 1:numel(opts.STIM_ELECTRODE)
            arrayMapIdx = find(arrayMap.chan == opts.STIM_ELECTRODE(stElec));
            if(~isempty(arrayMapIdx))
                rowElec = 11 - arrayMap.row(arrayMapIdx);
                colElec = arrayMap.col(arrayMapIdx);
                if(strcmpi(opts.STIM_ELECTRODE_LABEL,'string'))
                    if(numel(opts.STIM_ELECTRODE_COLOR) > 1)
                        text(colElec+0.25,rowElec+0.5,'S','FontWeight','bold','fontsize',20,'color',opts.STIM_ELECTRODE_COLOR{stElec});
                    else
                        text(colElec+0.25,rowElec+0.5,'S','FontWeight','bold','fontsize',20,'color',opts.STIM_ELECTRODE_COLOR);
                    end
                elseif(strcmpi(opts.STIM_ELECTRODE_LABEL,'box'))
                    if(numel(opts.STIM_ELECTRODE_COLOR) > 1)
                        rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',opts.STIM_ELECTRODE_COLOR{stElec},...
                            'FaceColor',opts.STIM_ELECTRODE_COLOR,'linewidth',0.1);
                    else
                        rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',opts.STIM_ELECTRODE_COLOR,...
                            'FaceColor',opts.STIM_ELECTRODE_COLOR,'linewidth',0.1);
                    end
                end
            end
        end
    end
    %% label recording electrode
    if(~isempty(opts.RECORDING_ELECTRODE))
        for recElec = 1:numel(opts.RECORDING_ELECTRODE)
            arrayMapIdx = find(arrayMap.chan == opts.RECORDING_ELECTRODE(recElec));
            if(~isempty(arrayMapIdx))
                rowElec = 11 - arrayMap.row(arrayMapIdx);
                colElec = arrayMap.col(arrayMapIdx);
                if(strcmpi(opts.RECORDING_ELECTRODE_LABEL,'string'))
                    if(numel(opts.RECORDING_ELECTRODE_COLOR)>1)
                        text(colElec+0.25,rowElec+0.5,'R','FontWeight','bold','fontsize',20,'color',opts.RECORDING_ELECTRODE_COLOR{recElec});
                    else
                        text(colElec+0.25,rowElec+0.5,'R','FontWeight','bold','fontsize',20,'color',opts.RECORDING_ELECTRODE_COLOR);
                    end
                elseif(strcmpi(opts.RECORDING_ELECTRODE_LABEL,'box'))
                    rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',opts.RECORDING_ELECTRODE_COLOR,...
                        'FaceColor',opts.RECORDING_ELECTRODE_COLOR,'linewidth',0.1);
                    if(numel(opts.RECORDING_ELECTRODE_COLOR)>1)
                        rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',opts.RECORDING_ELECTRODE_COLOR,...
                            'FaceColor',opts.RECORDING_ELECTRODE_COLOR{recElec},'linewidth',0.1);
                    else
                        rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',opts.RECORDING_ELECTRODE_COLOR,...
                            'FaceColor',opts.RECORDING_ELECTRODE_COLOR,'linewidth',0.1);
                    end
                end
            end
        end
    end

    ax = gca;
    set(ax,'Visible','off')
    if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_PREFIX,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
        saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'_arrayMap'));
    end

end



function [opts] = configureOpts(optsInput)

    opts.NUM_ROWS = 10;
    opts.NUM_COLS = 10;
    opts.STIM_ELECTRODE = [];
    opts.STIM_ELECTRODE_COLOR = {'k','r','b',[0,0.5,0],'m'};
    opts.STIM_ELECTRODE_LABEL = 'string';
    opts.RECORDING_ELECTRODE = [];
    opts.RECORDING_ELECTRODE_COLOR = 'k';
    opts.RECORDING_ELECTRODE_LABEL = 'string';
        
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.COLOR_MAP = 'white';
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
    

end % end function
