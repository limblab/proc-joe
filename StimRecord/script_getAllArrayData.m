%% combine array data's that exist within a folder structure

    single_pulse = 0;
    include_non_stim_data = 1;
    
    if(~single_pulse && include_non_stim_data)
        highest_folderpath{1} = 'E:\Data\Joseph\long_trains_array_data\amp\';
        search_word = 'arrayData';
    elseif(single_pulse)
        highest_folderpath{1} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Duncan\SinglePulse\';
        highest_folderpath{2} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\SinglePulse\';
        search_word = 'array_data';
    else 
        highest_folderpath{1} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Duncan\DblPulse_trains\';
        highest_folderpath{2} = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\DblPulse_trains\all_processed\';
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
            
            
            % get stim chan from filename
            chan_idx = strfind(files(file_num).name,'chan');
            underscore_idx = strfind(files(file_num).name,'_');
            chan_underscore_idx = underscore_idx(find(underscore_idx > chan_idx,1,'first'));
            
            stim_chan = num2str(files(file_num).name(chan_idx+4:chan_underscore_idx-1));
            % put arrayData into array_data_all;
            if(single_pulse || include_non_stim_data)
                num_units = numel(arrayData);
            else
                num_units = numel(array_data);
            end
            
            for arrayDataIdx = 1:num_units
                if(~single_pulse && include_non_stim_data)
                    % split out stim and non stim channels
                    monkey_underscore_idx = underscore_idx(1);
                    
                    if(arrayData{arrayDataIdx}.CHAN_REC == arrayData{arrayDataIdx}.CHAN_SENT{1})
                        array_data_stim_all{end+1} = arrayData{arrayDataIdx};
                        array_data_stim_all{end}.monkey = files(file_num).name(1:monkey_underscore_idx-1);
                    else
                        array_data_nonstim_all{end+1} = arrayData{arrayDataIdx};
                        chan_rec_list(end+1) = arrayData{arrayDataIdx}.CHAN_REC;
                        array_data_nonstim_all{end}.monkey = files(file_num).name(1:monkey_underscore_idx-1);
                    end
                    
                    
                elseif(single_pulse)
                    array_data_all{end+1} = arrayData{arrayDataIdx};
                    array_data_all{end}.CHAN_LIST = stim_chan;
                    array_data_all{end}.monkey = monkey_name{folder_num};
                else
                    array_data_all{end+1} = array_data{arrayDataIdx};
                    array_data_all{end}.stimData = array_data_all{end}.trial_num;
                    array_data_all{end}.numStims = array_data_all{end}.num_stims;
                    input_data_all{end+1} = input_data;
                    array_data_all{end}.CHAN_LIST = stim_chan;
                    array_data_all{end}.monkey = monkey_name{folder_num};
                end                
                
                
            end
        end
    end

    
    if(single_pulse)
        arrayData = array_data_all; % rename to match other scripts
    else
        array_data = array_data_all; % match with train script 
    end
   