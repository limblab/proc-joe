function [output_data] = predictRT(td,opts)

    %% configure opts and set default values
    opts = configureOpts(opts);
    output_data = [];
    
    num_bins = opts.WINDOW(2) - opts.WINDOW(1) + 1;
    opts.WINDOW(2) = opts.WINDOW(2) - mod(num_bins,opts.BIN_SIZE);
    
    if(isempty(opts.SPIKE_LIST))
        opts.SPIKE_LIST = 1:size(td(1).LeftS1_spikes,2);
    end
    if(isempty(opts.TRAIN_IDX))
        opts.TRAIN_IDX = 1:numel(td);
    end
    %% get x and y
    x = [];
    y = [];
    for t = 1:numel(td)
        y(end+1,1) = td(t).bin_size*(td(t).idx_movement_on - td(t).idx_goCueTime);
        idx = td(t).idx_goCueTime;
        spike_count = td(t).LeftS1_spikes(idx+opts.WINDOW(1):idx+opts.WINDOW(2),opts.SPIKE_LIST);
        % upbin if necessary
        for b = 1:opts.BIN_SIZE
            if(b == 1)
                binned_spike_count = spike_count(b:opts.BIN_SIZE:end,:);
            else
                binned_spike_count = binned_spike_count + spike_count(b:opts.BIN_SIZE:end,:);
            end
        end
        % each row of x is a constant, then all bins related to a unit, followed by all bins
        % related to the next unit
        x(end+1,:) = [1,reshape(binned_spike_count,1,numel(binned_spike_count))];
    end
    
    %% linear regression
    b = x(opts.TRAIN_IDX,:)\y(opts.TRAIN_IDX,:);
    
    %% configure output_data
    output_data.y_true = y;
%     output_data.y_train_idx = train_idx;
    output_data.y_pred = x*b;
    output_data.b = b;
end

function [opts] = configureOpts(optsInput)

   opts.SPIKE_LIST = [];
   opts.WINDOW = [0,9];
   opts.BIN_SIZE = 2; % # bins to combine
   opts.TRAIN_IDX = [];

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