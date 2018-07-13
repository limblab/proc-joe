function [output_data] = plotSpikesRT(td,opts)

    %% configure opts
    opts = configureOpts(opts,td);
    
    if(isempty(opts.SPIKE_LIST))
        opts.SPIKE_LIST = 1:size(td(1).LeftS1_spikes,2);
    end

    %% 
    rt = zeros(numel(td),1);
    num_spikes = zeros(numel(td),1);
    code = zeros(numel(td),1);
    
    exp_window = exp(-(mode([td.bin_size])*(0:opts.OFFSET))/opts.TAU);
    if(1 || isempty(opts.CORR_DIR))
        opts.CORR_DIR = ones(1,numel(opts.SPIKE_LIST));
    end
    
    baseline = [];
    for t = 1:numel(td)
        if t == 1
            baseline = (td(t).LeftS1_spikes(td(t).idx_goCueTime-30:td(t).idx_goCueTime-1,opts.SPIKE_LIST));
        else
            baseline = baseline + (td(t).LeftS1_spikes(td(t).idx_goCueTime-30:td(t).idx_goCueTime-1,opts.SPIKE_LIST));
        end
    end
    baseline = mean(baseline/numel(td));
    
    for t = 1:numel(td)
        rt(t) = td(t).bin_size*(td(t).idx_movement_on) - td(t).goCueTime;
        num_spikes(t) = sum(exp_window*(opts.CORR_DIR.*(td(t).LeftS1_spikes(td(t).idx_goCueTime:td(t).idx_goCueTime+opts.OFFSET,opts.SPIKE_LIST)-baseline)));
        if(opts.STIM)
            code(t) = td(t).stimCode;
        else
            code(t) = 10*td(t).bumpMagnitude;
        end
    end

    if(isempty(opts.STIM_CODES))
        opts.STIM_CODES = unique(code);
    end
    
    output_data.fits(numel(opts.STIM_CODES),1) = struct('stats',[],'fit_obj',[]);
    
    if(opts.MAKE_PLOTS)
        figure;
    end
    for s = 1:numel(opts.STIM_CODES)
        if(opts.MAKE_PLOTS)
            if(size(opts.COLORS,1) > 1)
    %             plot(rt(code == opts.STIM_CODES(s)),num_spikes(code == opts.STIM_CODES(s)),'.','markersize',12,'color',opts.COLORS(s,:));
                plot(num_spikes(isEqual(code',opts.STIM_CODES(s))),rt(isEqual(code',opts.STIM_CODES(s))),'.','markersize',12,'color',opts.COLORS(s,:));
            else
    %             plot(rt(code == opts.STIM_CODES(s)),num_spikes(code == opts.STIM_CODES(s)),'.','markersize',12,'color',opts.COLORS(1,:));
                plot(num_spikes(isEqual(code',opts.STIM_CODES(s))),rt(isEqual(code',opts.STIM_CODES(s))),'.','markersize',12,'color',opts.COLORS(1,:));
            end
        end
        if(sum(code == opts.STIM_CODES(s)) > opts.MIN_POINTS)
%             [output_data.fits(s).fit_obj,output_data.fits(s).stats] = fit(rt(code == opts.STIM_CODES(s)), num_spikes(code == opts.STIM_CODES(s)),'a*x+b');
            [output_data.fits(s).fit_obj,output_data.fits(s).stats] = fit(num_spikes(isEqual(code',opts.STIM_CODES(s))),rt(isEqual(code',opts.STIM_CODES(s))),'a*x+b');
 
%             output_data.corr(s) = corr(rt(code == opts.STIM_CODES(s)), num_spikes(code == opts.STIM_CODES(s)));
            output_data.corr(s) = corr(num_spikes(isEqual(code',opts.STIM_CODES(s))),rt(isEqual(code',opts.STIM_CODES(s))));
        end
        
        hold on
    end

    output_data.fit_all_codes = [];
%     [output_data.fit_all_codes.fit_obj,output_data.fit_all_codes.stats] = fit(rt,num_spikes,'a*x+b');
    [output_data.fit_all_codes.fit_obj,output_data.fit_all_codes.stats] = fit(num_spikes,rt,'a*x+b');
%     output_data.corr_all_codes = corr(rt,num_spikes);
    output_data.corr_all_codes = corr(num_spikes,rt);
    
    output_data.rt = rt;
    output_data.num_spikes = num_spikes;
    output_data.code = code;
    output_data.stim = opts.STIM;
end

function [opts] = configureOpts(optsInput,td)

    opts = [];
    opts.OFFSET = 9;
    opts.SPIKE_LIST = [];
    opts.TAU = 0.12; % exponential decaying window, set to large # for no decay
    opts.COLORS = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
    opts.STIM_CODES = [];
    opts.MIN_POINTS = 5;
    opts.STIM = 1;
    opts.MAKE_PLOTS = 1;
    
    opts.CORR_DIR = [];
    
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