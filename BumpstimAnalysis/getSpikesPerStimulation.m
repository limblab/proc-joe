function [ spikesPerStimulation, figureHandles ] = getSpikesPerStimulation( arrayData,opts )

    %% configure options and set default values
    opts = configureOpts(opts);
    
    spikesPerStimulation = cell(size(arrayData{1,1}.numStims,1),size(arrayData{1,1}.numStims,2));
        
    for chan = 1:size(arrayData{1,1}.numStims,1)
        for wave = 1:size(arrayData{1,1}.numStims,2)
            %% collect all unit data into two arrays (one for time, one for stim num)
            spikeTimes = [];
            stimNums = [];
            
            for unit = 1:numel(arrayData)
                spikeTimes = [spikeTimes, arrayData{unit}.spikeTrialTimes{chan,wave}];
                stimNums = [stimNums, arrayData{unit}.stimData{chan,wave}];
            end
            
            %% remove spikes that are not in the window
            spikeMask = spikeTimes > opts.WINDOW(1) & spikeTimes < opts.WINDOW(2);
            
            %% count how many spikes for each stimulation
            spikesPerStimulation{chan,wave} = zeros(arrayData{1,1}.numStims(chan,wave),1);
            
            for stim = 1:arrayData{1,1}.numStims(chan,wave)
                spikesPerStimulation{chan,wave}(stim) = sum(stimNums(spikeMask) == stim);
            end
            
        end
    end


        %% plot distribution for each channel/wave
    if(opts.MAKE_DISTRIBUTION_PLOT)
        for chan = 1:size(arrayData{1,1}.numStims,1)
            for wave = 1:size(arrayData{1,1}.numStims,2)
                figureHandles{chan,wave} = figure();
                [bC,bE] = histcounts(spikesPerStimulation{chan,wave});
                bE = bE(1:end-1) + (bE(2)-bE(1))/2;
                if(opts.PLOT_PERCENTAGE_STIMULATIONS)
                    bC = bC/arrayData{1,1}.numStims(chan,wave);
                    ylabel('Percent of stimulations')
                else
                    ylabel('Number of stimulations')
                end
                
                if(opts.PLOT_PERCENTAGE_UNITS)
                    bE = bE/numel(arrayData);
                    xlabel('Percent of units activated')
                else
                    xlabel('Number of evoked spikes')
                end
                
                 bar(bE,bC,'FaceColor',opts.DISTRIBUTION_FACECOLOR);
                
                formatForLee(gcf)
            end
        end
    end
    
        %% look at probability distribution assuming independence
    if(opts.MAKE_PDF_PLOT)
        for chan = 1:size(arrayData{1,1}.numStims,1)
            for wave = 1:size(arrayData{1,1}.numStims,2)
                
                %% get unit probabilities into one array
                probUnit = zeros(1,numel(arrayData));
                for unit = 1:numel(arrayData)
                    probUnit(unit) = arrayData{unit}.singleProb.meanSpikesEvoked_all(chan,wave);
                end
                
                %% compute probability of all combinations for a number of evoked spikes
                probCount = zeros(1+opts.PDF_EXTEND+max(spikesPerStimulation{chan,wave}),1);
                probCount(1) = prod(1-probUnit); % do k = 0 by hand
                for k = 1:max(spikesPerStimulation{chan,wave})+opts.PDF_EXTEND
                    unitsFiring = nchoosek(1:numel(arrayData),k);
                    probCombos = repmat(1-probUnit,size(unitsFiring,1),1);
                    for n = 1:size(unitsFiring,2) % set probability to firing prob
                        [row,col] = find(repmat(1:numel(arrayData),size(unitsFiring,1),1)==unitsFiring(:,n));
                        
                        probCombos((col-1).*size(probCombos,1) + row) = 1-probCombos((col-1).*size(probCombos,1) + row);
                    end
                    probCount(k+1) = sum(prod(probCombos,2));
                end
                
                if(opts.PLOT_PDF_ON_SAME_PLOT && opts.MAKE_DISTRIBUTION_PLOT)
                    figure(figureHandles{chan,wave});
                    hold on
                else
                    figure();  
                end
                bE = 0:max(spikesPerStimulation{chan,wave})+opts.PDF_EXTEND;
                bC = probCount;

                ylabel('Number of stimulations')
                
                if(opts.PLOT_PERCENTAGE_UNITS)
                    bE = bE/numel(arrayData);
                    xlabel('Percent of units activated')
                else
                    xlabel('Number of evoked spikes')
                end
                
                width = 1;
                if(opts.PLOT_PDF_ON_SAME_PLOT && opts.MAKE_DISTRIBUTION_PLOT)
                    width = 0.4;
                end
                bar(bE,bC,width,'FaceColor',opts.PDF_FACECOLOR)
                
                if(opts.PLOT_PDF_ON_SAME_PLOT && opts.MAKE_DISTRIBUTION_PLOT)
                    l = legend('recorded data','independence model');
                    set(l,'box','off')
                end
                
                xlabel('Number of evoked spikes per stimulation')
                ylabel('Probability')
                
                formatForLee(gcf)
            end
        end
    end
    
end

function [opts] = configureOpts(optsInput)

    opts.MAKE_DISTRIBUTION_PLOT = 1;
    opts.MAKE_PDF_PLOT = 1;
    opts.PLOT_PERCENTAGE_STIMULATIONS = 1;
    opts.PLOT_PERCENTAGE_UNITS = 1;
    opts.WINDOW = [1,5]/1000;
    
    opts.PDF_FACECOLOR = 'r';
    opts.DISTRIBUTION_FACECOLOR = 'k';
    opts.PLOT_PDF_ON_SAME_PLOT = 1;
    opts.PDF_EXTEND = 2;
    
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