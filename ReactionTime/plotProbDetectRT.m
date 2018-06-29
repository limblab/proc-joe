function [output_data] = plotProbDetectRT(data,opts)

    %% configure opts
    opts = configureOpts(opts);
    
    
    % for each dataset plot rt vs detect prob
    output_data.plots = figure;
    hold on
    for d = 1:numel(data)
        % use the psychometric fit to get a prob detect for each stim param
        stim_params = data{d}.opts.STIM_PARAMS;
        fit_stim = data{d}.psychometric_fit_stim;
        stim_detect_prob = fit_stim(1) + fit_stim(2)*erf(fit_stim(3)*(stim_params-fit_stim(4)));
        zero_val = max(0,fit_stim(1) + fit_stim(2)*erf(fit_stim(3)*(0-fit_stim(4))));
        stim_detect_prob = (stim_detect_prob - zero_val)/(1-zero_val)/(fit_stim(1)+fit_stim(2));
        % plot rt vs. stim_detect_prob
        for cueIdx = 1:numel(data{d}.cueInfo)
            if(data{d}.cueInfo(cueIdx).stimCode ~= -1 && data{d}.cueInfo(cueIdx).bumpMag == 0 && ~isempty(data{d}.cueInfo(cueIdx).rt))
                stimCode = data{d}.cueInfo(cueIdx).stimCode + 1;
                plot(stim_detect_prob(stimCode),mean(data{d}.cueInfo(cueIdx).rt),'.','color',opts.COLORS{d},'markersize',opts.MARKER_SIZE);
                hold on
            end
        end
    end
    xlim([0,1.0])
    
    % for each dataset, plot rt vs axis for the first dataset
    figure()
    for d = 1:numel(data)
        if(d == 1)
            base_fit = data{d}.psychometric_fit_stim;
            x_data = data{d}.opts.STIM_PARAMS;
        else
            % use the psychometric fit to get a prob detect for each stim param
            stim_params = data{d}.opts.STIM_PARAMS;
            fit_stim = data{d}.psychometric_fit_stim;
            stim_detect_prob = fit_stim(1) + fit_stim(2)*erf(fit_stim(3)*(stim_params-fit_stim(4)));
            zero_val = max(0,fit_stim(1) + fit_stim(2)*erf(fit_stim(3)*(0-fit_stim(4))));
            max_val = min(1,fit_stim(1) + fit_stim(2));
%             stim_detect_prob = (stim_detect_prob - zero_val)/(1-zero_val)/(fit_stim(1)+fit_stim(2));
            stim_detect_prob = (stim_detect_prob - zero_val)/(max_val-zero_val);
            x_data = erfinv((stim_detect_prob - base_fit(1))/base_fit(2))/base_fit(3) + base_fit(4);
        end
        % plot rt vs. stim_detect_prob
        for cueIdx = 1:numel(data{d}.cueInfo)
            if(data{d}.cueInfo(cueIdx).stimCode ~= -1 && data{d}.cueInfo(cueIdx).bumpMag == 0 && ~isempty(data{d}.cueInfo(cueIdx).rt))
                stimCode = data{d}.cueInfo(cueIdx).stimCode + 1;
                plot(x_data(stimCode),mean(data{d}.cueInfo(cueIdx).rt),'.','color',opts.COLORS{d},'markersize',opts.MARKER_SIZE);
                hold on
            end
        end
    end
    
    
    
end

function [opts] = configureOpts(optsInput)

    opts = [];
    opts.MARKER_SIZE =12;
    opts.COLORS = {'r',[0 0.5 0],'b'}'
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
    

end