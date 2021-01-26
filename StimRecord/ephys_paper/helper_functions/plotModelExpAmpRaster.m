function [output_data] = plotModelExpAmpRaster(data, input_data)

    
    % send in x_data, y_data, optsPlot, optsSave
    optsPlot = [];
    optsPlot.DIVIDING_LINES = []; % list of numbers -- separate amplitude trials
    optsPlot.DIVIDING_LINES_COLORS = {};
    optsPlot.DUKE_BLANKED_REGION = []; % [min, max] for each region in dividing lines
    optsPlot.BLACK_BLANKED_REGION = []; % [min, max] for each region in dividing lines
    optsSave = [];
    x_data = [];
    y_data = [];
    
    % get x_data and y_data -- do for each amplitude
    trial_counter = 1;
    if(input_data.plot_amp)
        for i_amp = 1:numel(input_data.amp_list)
            amp_idx = find(input_data.amp_list(i_amp) == [data.STIM_PARAMETERS.amp1],1,'first');

            for i_stim = 1:data.numStims(amp_idx)
                mask = data.stimData{amp_idx} == i_stim;
                if(sum(mask) > 0)
                    x_data(end+1:end+sum(mask)) =  data.spikeTrialTimes{amp_idx}(mask)*1000; % convert to ms
                    y_data(end+1:end+sum(mask)) = trial_counter*ones(1,sum(mask));
                end
                trial_counter = trial_counter + 1;
            end
            optsPlot.DIVIDING_LINES(end+1) = trial_counter - 0.5;
            optsPlot.DIVIDING_LINES_COLORS{end+1} = getColorFromList(1,1);
            
            % get blanked region
            if(isfield(input_data,'duke_blank_times'))
                optsPlot.DUKE_BLANKED_REGION(end+1,:) = input_data.duke_blank_times(amp_idx,:);
            end
            
            % get blanked region
            if(isfield(input_data,'black_blank_times'))
                optsPlot.BLACK_BLANKED_REGION(end+1,:) = input_data.black_blank_times(amp_idx,:);
            end
        end
        optsPlot.DIVIDING_LINES = optsPlot.DIVIDING_LINES(1:end-1); % remove last entry
    else % plot polarity
            % check to see if 2 50uA entries exist, and if they have
            % different polarities
        pol_idx = [find([data.STIM_PARAMETERS.amp1] == 50 & [data.STIM_PARAMETERS.polarity]==0),...
            find([data.STIM_PARAMETERS.amp1] == 50 & [data.STIM_PARAMETERS.polarity]==1)];
        if(numel(pol_idx) == 2)
            for i_pol = 1:numel(pol_idx) % cathodic is plotted first
                for i_stim = 1:data.numStims(pol_idx(i_pol))
                    mask = data.stimData{pol_idx(i_pol)} == i_stim;
                    if(sum(mask) > 0)
                        x_data(end+1:end+sum(mask)) = data.spikeTrialTimes{pol_idx(i_pol)}(mask)*1000;
                        y_data(end+1:end+sum(mask)) = trial_counter+ones(1,sum(mask));
                    end
                    trial_counter = trial_counter + 1;
                end
                optsPlot.DIVIDING_LINES(end+1) = trial_counter - 0.5;
                optsPlot.DIVIDING_LINES_COLORS{end+1} = getColorFromList(1,1);
                
            end
            optsPlot.DIVIDING_LINES = optsPlot.DIVIDING_LINES(1:end-1); % remove last entry
        end
    end
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
    
    
    if(~isempty(x_data))
        output_data = plotRasterLIB(x_data,y_data,optsPlot,optsSave);
    end

end