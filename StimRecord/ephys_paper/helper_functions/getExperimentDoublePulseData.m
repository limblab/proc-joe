function [output_data] = getExperimentDoublePulseData(input_data)
% input_data contains
    % home_computer : bool, which sets folderpaths
    
% output_data contains 
    % array_data : list of array_data structs for each neuron
    
    
    if(input_data.home_computer)
        highest_folderpath{1} = 'D:\Lab\Data\stim_ephys_paper\dbl_pulse_short_trains_data\Duncan\';
        highest_folderpath{2} = 'D:\Lab\Data\stim_ephys_paper\dbl_pulse_short_trains_data\Han\';
        search_word = 'arrayData';
    else
        error('file path not implemented for laptop yet');
        search_word = 'arrayData';
    end
    
    monkey_name{1} = 'Duncan';
    monkey_name{2} = 'Han';

    array_data_all = {};
    array_data_stim_all = {};
    array_data_nonstim_all = {};
    input_data_all = {};
    chan_rec_list = [];
    
    pwd = cd;
    for folder_num = 1:numel(highest_folderpath)
        cd(highest_folderpath{folder_num})
        files = dir([highest_folderpath{folder_num} '**\*',search_word,'*']);

        for file_num = 1:numel(files)
            % load files(f)
            load([files(file_num).folder,'\' files(file_num).name]);
            
            underscore_idx = strfind(files(file_num).name,'_');

            % put arrayData into array_data_all;
            num_units = numel(array_data);
            
            for arrayDataIdx = 1:num_units
                % get stim chan from filename
                chan_idx = strfind(files(file_num).name,'chan');
                chan_underscore_idx = underscore_idx(find(underscore_idx > chan_idx,1,'first'));

                stim_chan = num2str(files(file_num).name(chan_idx+4:chan_underscore_idx-1));

                array_data_all{end+1} = array_data{arrayDataIdx};
                array_data_all{end}.CHAN_LIST = stim_chan;
                array_data_all{end}.monkey = monkey_name{folder_num};
            end
        end
    end    
    
    % package outputs
    output_data.array_data = array_data_all;

end