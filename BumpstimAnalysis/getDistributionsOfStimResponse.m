function [arrayDataFits] = getDistributionsOfStimResponse( arrayData,opts )
% 

    opts = configureOpts(opts);

    %% get data for each window
    means = zeros(size(opts.WINDOW,1),size(arrayData,2),size(arrayData{1}.bC,1),size(arrayData{1}.bC,2));
    vars = zeros(size(opts.WINDOW,1),size(arrayData,2),size(arrayData{1}.bC,1),size(arrayData{1}.bC,2));
    
    for windowIdx = 1:size(opts.WINDOW,1) % for each window
        for arrayIdx = 1:size(arrayData,2) % for each unit
            for chan = 1:size(arrayData{arrayIdx}.bC,1) % for each channel stimulated
                for waveNum = 1:size(arrayData{arrayIdx}.bC,2) % for each waveform stimulated with
                    binIdx = [find(abs(arrayData{arrayIdx}.bE{chan,waveNum}-opts.WINDOW(windowIdx,1))<1E-4), find(abs(arrayData{arrayIdx}.bE{chan,waveNum}-opts.WINDOW(windowIdx,2))<1E-4)];
                    means(windowIdx,arrayIdx,chan,waveNum) = mean(arrayData{arrayIdx}.bC{chan,waveNum}(binIdx(1):binIdx(2)));
                    vars(windowIdx,arrayIdx,chan,waveNum) = var(arrayData{arrayIdx}.bC{chan,waveNum}(binIdx(1):binIdx(2)));
                end % end for wave
            end % end for chan
        end % end for unit
    end % end for window
    
    %% plot data for each window with a best fit line going through (0,0)
    windowIdxs = 1:1:size(opts.WINDOW,1);
    if(~isempty(opts.BASELINE_WINDOW_IDX))
        windowIdxs(opts.BASELINE_WINDOW_IDX) = 1;
        windowIdxs(1) = opts.BASELINE_WINDOW_IDX;
    end
    for windowIdx = windowIdxs
        % get data and compress across all units, channels, waveforms sent
        meanWindow = reshape(means(windowIdx,:,:,:),numel(means(windowIdx,:,:,:)),1);
        varWindow = reshape(vars(windowIdx,:,:,:),numel(vars(windowIdx,:,:,:)),1);
        
        % plot variance vs mean
        
        [F{windowIdx},gof{windowIdx}] = fit(meanWindow,varWindow,'a*x');
        xFit = [0,max(meanWindow)];
        if(isempty(opts.BASELINE_WINDOW_IDX))
            yFit = F{windowIdx}.a*xFit;
        else
            yFit = F{opts.BASELINE_WINDOW_IDX}.a*xFit;
        end
        
        if(~isempty(opts.BASELINE_WINDOW_IDX) && opts.PROJECT_TO_SLOPE_1)
            xFit = xFit*F{opts.BASELINE_WINDOW_IDX}.a;
            meanWindow = meanWindow*F{opts.BASELINE_WINDOW_IDX}.a;
        end
        figure();
        if(opts.COLOR_MARKERS_RESPONSE && opts.BASELINE_WINDOW_IDX)
            meanBaseline = reshape(means(opts.BASELINE_WINDOW_IDX,:,:,:),numel(means(opts.BASELINE_WINDOW_IDX,:,:,:)),1);
            if(opts.PROJECT_TO_SLOPE_1)
                meanBaseline = meanBaseline*F{opts.BASELINE_WINDOW_IDX}.a;
            end
            ratios = meanWindow./meanBaseline;
            
            for r = 1:numel(ratios)
                plot(meanWindow(r),varWindow(r),'.','markersize',opts.MARKER_SIZE,'color',opts.COLORS(floor(size(opts.COLORS,1)*ratios(r)/max(ratios)),:));
                hold on
            end
        else
            plot(meanWindow,varWindow,'.','markersize',opts.MARKER_SIZE);
        end
        hold on
        plot(xFit,yFit,'k')
        
    end
    
    arrayDataFits.F = F;
    arrayDataFits.gof = gof;
end

function [opts] = configureOpts(optsInput)

    opts.WINDOW = [];
    opts.MARKER_SIZE = 12;
    opts.BASELINE_WINDOW_IDX = [];
    opts.PROJECT_TO_SLOPE_1 = 0;
    opts.COLOR_MARKERS_RESPONSE = 0;
    opts.COLORS = colormap(jet);

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