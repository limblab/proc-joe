function [amp_freq_data, intermittent_data] = getExperimentLongTrainData(input_data)

    if(~input_data.home_computer)
%         highest_folderpath{1} = 'D:\Lab\Data\StimArtifact\long_trains_array_data\';
        erorr('folderpath for laptop not yet specified');
    else
        highest_folderpath = 'D:\Lab\Data\StimArtifact\long_trains_array_data\';
    end
    
    search_word = 'arrayData';
    
    monkey_name{1} = 'Duncan';
    monkey_name{2} = 'Han';

    amp_freq_stim = {};
    amp_freq_nonstim = {};
    amp_freq_chan_rec = [];
    
    inter_180_stim = {};
    inter_180_nonstim = {};
    inter_180_chan_rec = [];
    
    inter_130_stim = {};
    inter_130_nonstim = {};
    inter_130_chan_rec = [];
    
    pwd = cd;
    % load in amp_freq data
    files = dir([highest_folderpath 'amp_freq\pulse_times\*',search_word,'*']);

    for file_num = 1:numel(files)
        % load files(f)
        load([files(file_num).folder,'\' files(file_num).name]);

        % put arrayData into array_data_all;
        [stim_all_temp,nonstim_all_temp,chan_rec_temp] = loadInData(arrayData,files(file_num).name);

        amp_freq_stim(end+1:end+numel(stim_all_temp)) = stim_all_temp;
        amp_freq_nonstim(end+1:end+numel(nonstim_all_temp)) = nonstim_all_temp;
        amp_freq_chan_rec(end+1:end+numel(nonstim_all_temp)) = chan_rec_temp;
    end

    % load in intermittent data (180 Hz)
    files = dir([highest_folderpath 'intermittent_180Hz\pulse_times\*',search_word,'*']);

    for file_num = 1:numel(files)
        % load files(f)
        load([files(file_num).folder,'\' files(file_num).name]);

        % put arrayData into array_data_all;
        [stim_all_temp,nonstim_all_temp,chan_rec_temp] = loadInData(arrayData,files(file_num).name);

        inter_180_stim(end+1:end+numel(stim_all_temp)) = stim_all_temp;
        inter_180_nonstim(end+1:end+numel(nonstim_all_temp)) = nonstim_all_temp;
        inter_180_chan_rec(end+1:end+numel(nonstim_all_temp)) = chan_rec_temp;
    end
    
    % load in intermittent data (130 Hz)
    files = dir([highest_folderpath 'intermittent_130Hz\pulse_times\*',search_word,'*']);
        
    for file_num = 1:numel(files)
        % load files(f)
        load([files(file_num).folder,'\' files(file_num).name]);

        % put arrayData into array_data_all;
        [stim_all_temp,nonstim_all_temp,chan_rec_temp] = loadInData(arrayData,files(file_num).name);

        inter_130_stim(end+1:end+numel(stim_all_temp)) = stim_all_temp;
        inter_130_nonstim(end+1:end+numel(nonstim_all_temp)) = nonstim_all_temp;
        inter_130_chan_rec(end+1:end+numel(nonstim_all_temp)) = chan_rec_temp;
    end
    
    % package outputs
    amp_freq_data.stim_chan = amp_freq_stim;
    amp_freq_data.nonstim_chan = amp_freq_nonstim;
    amp_freq_data.nonstim_chan_rec = amp_freq_chan_rec;
    
    intermittent_data.high_freq.stim_chan = inter_180_stim;
    intermittent_data.high_freq.nonstim_chan = inter_180_nonstim;
    intermittent_data.high_freq.nonstim_chan_rec = inter_180_chan_rec;
    
    intermittent_data.low_freq.stim_chan = inter_130_stim;
    intermittent_data.low_freq.nonstim_chan = inter_130_nonstim;
    intermittent_data.low_freq.nonstim_chan_rec = inter_130_chan_rec;
end




% function to load in data and separate into stim_chan and non_stim_chan
% results
function [stim_all, nonstim_all,chan_rec] = loadInData(arrayData,filename)
    stim_all = {}; nonstim_all = {}; chan_rec = [];          
    
    underscore_idx = strfind(filename,'_');

    for arrayDataIdx = 1:numel(arrayData)
        % split out stim and non stim channels
        monkey_underscore_idx = underscore_idx(1);

        if(arrayData{arrayDataIdx}.CHAN_REC == arrayData{arrayDataIdx}.CHAN_SENT{1})
            stim_all{end+1} = arrayData{arrayDataIdx};
            stim_all{end}.monkey = filename(1:monkey_underscore_idx-1);
        else
            nonstim_all{end+1} = arrayData{arrayDataIdx};
            chan_rec(end+1) = arrayData{arrayDataIdx}.CHAN_REC;
            nonstim_all{end}.monkey = filename(1:monkey_underscore_idx-1);
        end


    end


end
