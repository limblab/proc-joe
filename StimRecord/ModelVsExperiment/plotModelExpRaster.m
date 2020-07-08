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
    
    for i_amp = 1:numel(data.spikeTrialTimes)
        for i_stim = 1:data.numStims(i_amp)
            mask = data.stimData{i_amp} == i_stim;
            if(sum(mask > 0))
                x_data(end+1:end+sum(mask)) =  data.spikeTrialTimes{i_amp}(mask)*1000; % convert to ms
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
    optsPlot.MARKER_STYLE = 'line';
    % populate optsSave
    optsSave = [];
    

    output_data = plotRasterLIB(x_data,y_data,optsPlot,optsSave);


end