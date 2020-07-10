function [output_data] = plotModelExpRaster(data, input_data)

    
    % send in x_data, y_data, optsPlot, optsSave
    optsPlot = [];
    optsPlot.DIVIDING_LINES = []; % list of numbers -- separate amplitude trials
    optsPlot.DIVIDING_LINES_COLORS = {};
    optsSave = [];
    x_data = [];
    y_data = [];
    
    % get x_data and y_data -- do for each amplitude
    trial_counter = 1;
    
    for i_amp = 1:numel(input_data.amp_list)
        amp_idx = find(input_data.amp_list(i_amp) == [data.STIM_PARAMETERS.amp1],1,'first');
        
        for i_stim = 1:data.numStims(amp_idx)
            mask = data.stimData{amp_idx} == i_stim;
            if(sum(mask > 0))
                x_data(end+1:end+sum(mask)) =  data.spikeTrialTimes{amp_idx}(mask)*1000; % convert to ms
                y_data(end+1:end+sum(mask)) = trial_counter*ones(1,sum(mask));
            end
            trial_counter = trial_counter + 1;
        end
        optsPlot.DIVIDING_LINES(end+1) = trial_counter - 0.5;
        optsPlot.DIVIDING_LINES_COLORS{end+1} = getColorFromList(1,1);
    end
    optsPlot.DIVIDING_LINES = optsPlot.DIVIDING_LINES(1:end-1); % remove last entry
    
    % populate optsPlot
    optsPlot.X_LABEL = 'Time after stim offset (ms)';
    optsPlot.Y_LABEL = 'Stimuli';
    optsPlot.PLOT_ZERO_LINE = 0;
    optsPlot.X_LIMITS = input_data.x_lim;
    optsPlot.MARKER_STYLE = input_data.marker_style;
    if(input_data.is_model)
        optsPlot.LINE_WIDTH = 2;
    else
        optsPlot.LINE_WIDTH = 3;
    end
    % populate optsSave
    optsSave = [];
    
    
    
    output_data = plotRasterLIB(x_data,y_data,optsPlot,optsSave);


end