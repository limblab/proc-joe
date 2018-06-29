function [plots] = plotRasterPSTHRT(td,opts)

    %% configure opts
    opts = configureOpts(opts);

    %% plot bump trials -- sort rasters into groups based on bump mag, then sort within each group based on rt
    if(opts.PLOT_BUMP)
        plot_order = 1:1:numel(td); % idx in td
        plot_order(~[td(plot_order).isBumpTrial]) = []; % remove non bump trials
        plot_order([td(plot_order).bumpMagnitude]==0) = []; % remove 0 mag bump trials (sanity check)
        plot_order(isnan([td(plot_order).idx_movement_on])) = []; % remove nan movement on trials
        [~,bump_sort] = sort([td(plot_order).bumpMagnitude]);% sort by bump mag
        plot_order = plot_order(bump_sort);

        % within each unique bump mag, sort by rt
        bump_mags = unique([td(plot_order).bumpMagnitude]);
        start_idx = 1;
        dividing_lines = [];
        for bm = 1:numel(bump_mags)
            sub_plot_order = plot_order(find(isEqual([td(plot_order).bumpMagnitude],bump_mags(bm))));
            [~,sub_plot_sort] = sort([td(sub_plot_order).idx_movement_on]-[td(sub_plot_order).idx_goCueTime]);
            plot_order(start_idx:start_idx+numel(sub_plot_order)-1) = sub_plot_order(sub_plot_sort);
            start_idx = start_idx+numel(sub_plot_order);

            dividing_lines = [dividing_lines,start_idx-0.5];
        end

        % plot the raster based on plot order, separate bumps with a line
        if(opts.PLOT_RASTER)
            for unit = 1:size(td(1).LeftS1_spikes,2)
                f = figure();
                xData_raster = [];
                yData_raster = [];
                for trial = 1:numel(plot_order)
                    trial_idx = plot_order(trial);
                    if(~iempty(td(trial_idx).LeftS1_ts{unit}))
                        xData_raster = [xData_raster;td(trial_idx).LeftS1_ts{unit}-td(trial_idx).goCueTime+opts.EXTRA_TIME(1)];
                        yData_raster = [yData_raster;trial*ones(numel(td(trial_idx).LeftS1_ts{unit}),1)];
                    end
                end
                % deal with optsSave and optsPlot
                optsSave.FIGURE_SAVE = opts.SAVE_FIGURES;
                f.Name = strcat(opts.FIGURE_PREFIX,'_nn',num2str(unit),'_chan',num2str(td(1).LeftS1_unit_guide(1)),'_unit',num2str(td(1).LeftS1_unit_guide(2)),'_bumpRaster');
                optsSave.FIGURE_NAME = f.Name;
                optsSave.FIGURE_DIR = opts.FIGURE_PATH;

                optsPlot.DIVIDING_LINES = dividing_lines;
                optsPlot.X_LIMITS = opts.X_LIMITS;
                optsPlot.MAKE_FIGURE = 0;
                plotRasterLIB(xData_raster,yData_raster,optsPlot,optsSave);
            end
        end
        % plot the PSTH, a line for each bump mag
        dividing_marks = [1,dividing_lines + 0.5];
        offset = floor(opts.X_LIMITS/td(1).bin_size);
        if(opts.PLOT_PSTH)
            for unit = 1:size(td(1).LeftS1_spikes,2)
                num_trials = 0;
                f = figure();
                xData_PSTH = repmat(([offset(1):1:offset(2)]*td(1).bin_size)',1,numel(dividing_marks)-1);
                yData_PSTH = zeros(size(xData_PSTH));
                for dm = 1:numel(dividing_marks)-1
                    trials = plot_order(dividing_marks(dm):dividing_marks(dm+1)-1);
                    for t = 1:numel(trials)
                        if(~isnan(td(trials(t)).LeftS1_spikes(1,unit)))
                            yData_PSTH(:,dm) = yData_PSTH(:,dm) + td(trials(t)).LeftS1_spikes(td(trials(t)).idx_goCueTime+offset(1):td(trials(t)).idx_goCueTime+offset(2),unit);
                            num_trials = num_trials+1;
                        end
                    end
                    % normalize by number of trials
                    yData_PSTH(:,dm) = yData_PSTH(:,dm)/num_trials;
                end
                % deal with optsSave and optsPlot
                optsSave.FIGURE_SAVE = opts.SAVE_FIGURES;
                f.Name = strcat(opts.FIGURE_PREFIX,'_nn',num2str(unit),'_chan',num2str(td(1).LeftS1_unit_guide(1)),'_unit',num2str(td(1).LeftS1_unit_guide(2)),'_bumpPSTH');
                optsSave.FIGURE_NAME = f.Name;
                optsSave.FIGURE_DIR = opts.FIGURE_PATH;

                optsPlot.DIVIDING_LINES = dividing_lines;
                optsPlot.X_LIMITS = opts.X_LIMITS;
                optsPlot.MAKE_FIGURE = 0;
                optsPlot.NUM_PLOTS = size(yData_PSTH,2);
                optsPlot.BAR_STYLE = 'line';
                plotPSTHLIB(xData_PSTH,yData_PSTH,optsPlot,optsSave);
            end
        end
    end
    %% plot stim trials -- sort rasters into groups based on stim code, then sort within each group based on rt
    if(opts.PLOT_STIM)
        plot_order = 1:1:numel(td); % idx in td
        plot_order(~[td(plot_order).isStimTrial]) = []; % remove non stim trials
        plot_order([td(plot_order).stimCode]==-1) = []; % remove bad stim trials (sanity check)
        plot_order(isnan([td(plot_order).idx_movement_on])) = []; % remove nan movement on trials
        
        [~,stim_sort] = sort([td(plot_order).stimCode]);% sort by stim code
        plot_order = plot_order(stim_sort);

        % within each unique stim code, sort by rt
        stim_codes = unique([td(plot_order).stimCode]);
        start_idx = 1;
        dividing_lines = [];
        for sc = 1:numel(stim_codes)
            sub_plot_order = plot_order(find(isEqual([td(plot_order).stimCode],stim_codes(sc))));
            [~,sub_plot_sort] = sort([td(sub_plot_order).idx_movement_on]-[td(sub_plot_order).idx_goCueTime]);
            plot_order(start_idx:start_idx+numel(sub_plot_order)-1) = sub_plot_order(sub_plot_sort);
            start_idx = start_idx+numel(sub_plot_order);

            dividing_lines = [dividing_lines,start_idx-0.5];
        end

        % plot the raster based on plot order, separate bumps with a line
        if(opts.PLOT_RASTER)
            for unit = 1:size(td(1).LeftS1_spikes,2)
                f = figure();
                xData_raster = [];
                yData_raster = [];
                for trial = 1:numel(plot_order)
                    trial_idx = plot_order(trial);
                    xData_raster = [xData_raster;td(trial_idx).LeftS1_ts{unit}-td(trial_idx).goCueTime+opts.EXTRA_TIME(1)];
                    yData_raster = [yData_raster;trial*ones(numel(td(trial_idx).LeftS1_ts{unit}),1)];
                end
                % deal with optsSave and optsPlot
                optsSave.FIGURE_SAVE = opts.SAVE_FIGURES;
                f.Name = strcat(opts.FIGURE_PREFIX,'_nn',num2str(unit),'_chan',num2str(td(1).LeftS1_unit_guide(1)),'_unit',num2str(td(1).LeftS1_unit_guide(2)),'_stimRaster');
                optsSave.FIGURE_NAME = f.Name;
                optsSave.FIGURE_DIR = opts.FIGURE_PATH;

                optsPlot.DIVIDING_LINES = dividing_lines;
                optsPlot.X_LIMITS = opts.X_LIMITS;
                optsPlot.MAKE_FIGURE = 0;
                plotRasterLIB(xData_raster,yData_raster,optsPlot,optsSave);
            end
        end
        % plot the PSTH, a line for each stim code
        dividing_marks = [1,dividing_lines + 0.5];
        offset = floor(opts.X_LIMITS/td(1).bin_size);
        if(opts.PLOT_PSTH)
            for unit = 1:size(td(1).LeftS1_spikes,2)
                f = figure();
                xData_PSTH = repmat(([offset(1):1:offset(2)]*td(1).bin_size)',1,numel(dividing_marks)-1);
                yData_PSTH = zeros(size(xData_PSTH));
                for dm = 1:numel(dividing_marks)-1
                    trials = plot_order(dividing_marks(dm):dividing_marks(dm+1)-1);
                    for t = 1:numel(trials)
                        yData_PSTH(:,dm) = yData_PSTH(:,dm) + td(trials(t)).LeftS1_spikes(td(trials(t)).idx_goCueTime+offset(1):td(trials(t)).idx_goCueTime+offset(2),unit);
                    end
                    % normalize by number of trials
                    yData_PSTH(:,dm) = yData_PSTH(:,dm)/numel(trials);
                end
                % deal with optsSave and optsPlot
                optsSave.FIGURE_SAVE = opts.SAVE_FIGURES;
                f.Name = strcat(opts.FIGURE_PREFIX,'_nn',num2str(unit),'_chan',num2str(td(1).LeftS1_unit_guide(1)),'_unit',num2str(td(1).LeftS1_unit_guide(2)),'_stimPSTH');
                optsSave.FIGURE_NAME = f.Name;
                optsSave.FIGURE_DIR = opts.FIGURE_PATH;

                optsPlot.DIVIDING_LINES = dividing_lines;
                optsPlot.X_LIMITS = opts.X_LIMITS;
                optsPlot.MAKE_FIGURE = 0;
                optsPlot.NUM_PLOTS = size(yData_PSTH,2);
                optsPlot.BAR_STYLE = 'line';
                plotPSTHLIB(xData_PSTH,yData_PSTH,optsPlot,optsSave);
            end
        end    
    end
end


function [opts] = configureOpts(optsInput)

    opts = [];
    opts.X_LIMITS = [-1,1];
    opts.EXTRA_TIME = [0.2,0.2]; % default in parseFileByTrial
    
    opts.FIGURE_PREFIX = '';
    opts.FIGURE_PATH = '';
    opts.SAVE_FIGURES = 1;
    
    opts.PLOT_RASTER = 1;
    opts.PLOT_PSTH = 1;
    
    opts.PLOT_STIM = 1;
    opts.PLOT_BUMP = 1;
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldNames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldNames)
           if(isfield(opts,inputFieldNames{fn}))
               opts.(inputFieldNames{fn}) = optsInput.(inputFieldNames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    

end